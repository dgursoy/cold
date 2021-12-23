#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <ctype.h>


#ifndef MAX_FILE_LENGTH
#define MAX_FILE_LENGTH 2048
#endif

#define MAX_NUM_IMAGE_FILES 20			/* maximum number of image files allowed (max number of detectors too) */
#define EXIT_WITH_HELP { for(i=0;help[i][0];i++) fprintf(stderr,"%s\n",help[i]); exit(1); }

#define MIN(A,B)	( (A)<(B)?(A):(B) )
#define MAX(A,B)	( (A)>(B)?(A):(B) )
#define CHECK_FREE(A)   { if(A) free(A); (A)=NULL;}


int main (int argc, const char * argv[]) {
	/* after running peaksearch and pixels2qs on multiple images, combines
	 *	the PeaksFile's with Q's into one file for indexing.
	 */
	static char *help[] = {"combinePeaks file1 [file2 file3 ...] outputFile", "\tno switches", "" };
	long	i;

	char	inFiles[MAX_NUM_IMAGE_FILES][MAX_FILE_LENGTH];	/* holds names of the input files */
	char	outFile[MAX_FILE_LENGTH];						/* holds name of the output file */
	char	*out=NULL;										/* buffer to hold output for outFile */
	char	*in=NULL;										/* buffer to hold contents of one input file */
	char	*header=NULL;									/* holds the file header (from first input file) */
	int		Nin;											/* number of input files */
	size_t	inLen=0, outLen=0;								/* current allocated lengths of buffers */
	size_t	fsize;											/* length of current input file */
	size_t	hlen=0;											/* lenght of header */
	char	N_Ghat[1024];									/* the $N_Ghat line */
	FILE	*f;
	char	*p;
	char	*p0,*p1;										/* points to start and end of data block */
	int		N;												/* number of peaks in one file */
	int		Npeaks=0;										/* total number of peaks in all files */
	int		first=1;

	if (argc<3) { fprintf(stderr,"ERROR: not enough inputs\n"); EXIT_WITH_HELP }
	for (i=1;i<(argc-1);i++) {
		if (argv[i][0]=='-') { fprintf(stderr,"ERROR: no switches allowed\n"); EXIT_WITH_HELP }
		strncpy(inFiles[i-1],argv[i],MAX_FILE_LENGTH);		/* an input file */
	}
	Nin = i-1;											/* number of input files */
	strncpy(outFile,argv[argc-1],MAX_FILE_LENGTH);		/* last arg is output file */

	for (i=0;i<Nin;i++) {								/* accumulate the input information */
/*		if ( !(f=fopen(inFiles[i],"r")) ) {fprintf(stderr,"ERROR: cannot open file '%s' to read\n",inFiles[i]); exit(1);} */
		if ( !(f=fopen(inFiles[i],"r")) ) continue;		/* This image had no peaks, look at the output from the other images */
		fseek(f, 0, SEEK_END);							/* seek to end of file, get size and reset to start */
		fsize = ftell(f);
		fseek(f, 0, SEEK_SET);
		if (fsize<1) continue;
		inLen = MAX(fsize+10,inLen);
		in = realloc(in,inLen);							/* allocate or extend input buffer */
		if (!in) {fprintf(stderr,"ERROR: cannot allocate %lu bytes for file '%s' to read\n",inLen,inFiles[i]); exit(1);}
		fread(in,sizeof(char),fsize,f);
		in[fsize] = '\0';
		fclose(f);
		if (!(p = strstr(in,"\n$N_Ghat+Intens")+1)) continue;
		if ((N = atoi(strlen("\n$N_Ghat+Intens")+p)) < 1) continue;

		if (!header) {									/* no header yet, read from here */
			hlen = p - in + 1;							/* length of header */
			header = calloc(hlen, sizeof(char)); 
			if (!header) {fprintf(stderr,"ERROR: cannot allocate %lu bytes for header\n",hlen); exit(1);}
			memcpy(header,in,hlen-1);
			header[hlen-1] = '\0';
		}

		if (!(p0=strchr(p,'\n'))) continue;				/* p0 is the start of the data block */
		p0++;

		p1 = in+fsize;									/* find p1, the end of the data block */
		for (p1=in+fsize; p1>in && !isdigit(*p1); p1--) {}
		*(++p1) = '\n';
		*(++p1) = '\0';

		outLen += strlen(p0);
		out = realloc(out,outLen+1);					/* allocate or extend output buffer */
		if (!in) {fprintf(stderr,"ERROR: cannot allocate %lu bytes for file '%s' to read\n",inLen,inFiles[i]); exit(1);}
		if (first) {out[0] = '\0'; first=0;}
		strcat(out,p0); 

		Npeaks += N;
	}
	CHECK_FREE(in)
	if (!out || !header) {fprintf(stderr,"ERROR: no peaks to write\n");exit(1);}
	sprintf(N_Ghat,"$N_Ghat+Intens \t%d\t\t// number of G^ vectors\n",Npeaks);


	if ( !(f=fopen(outFile,"w")) ) {fprintf(stderr,"ERROR: cannot open file '%s' to write results\n",outFile);exit(1);}
	fwrite(header,sizeof(char),hlen-1,f); 
	fwrite(N_Ghat,sizeof(char),strlen(N_Ghat),f); 
	fwrite(out,sizeof(char),strlen(out),f); 
	fclose(f);

	CHECK_FREE(header)
	CHECK_FREE(out)
	return 0;
}
