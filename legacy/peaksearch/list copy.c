/**********************************************************

	Data Structure and Routines for a doubly linked list.

/**********************************************************/

#include "list.h"

List* list_new(void){

	List* l = malloc(sizeof(List));
	if (!l) exit(ENOMEM);		/* Not enough space. */
	l->size = 0;
	l->head = EMPTY_NODE;
	l->tail = EMPTY_NODE;
	
	return l;

}

/*list_append takes a pointer to the data */
/*we include a pointer to the list so that we can update length*/
void list_append(List* l, void* value){

	/*create our new node*/
	ListNode* n = malloc(sizeof(ListNode));
	if (!n) exit(ENOMEM);		/* Not enough space. */
	/*set aside the required amount of space, and copy the data into it*/
	n->value = value;
	
	n->next = EMPTY_NODE;
	n->last = l->tail;
	
	/*if this list has no nodes in it yet...*/
	if(l->head == EMPTY_NODE) {
	
		/*set the list's head and tail to point to this node*/
		l->head = n;
		l->tail = n;
	
	} else {
	
		/*set the tail's 'next' property to point to this node, now, then set list's 'tail' property to point to n, too*/
		l->tail->next = n;
		l->tail = n;
		
	}
	
	/*increment the size of the list*/
	l->size++;
	

}

/*list_remove removes the given node from the list, and links the ones before and after it*/
/*we include a pointer to the list so that we can update length, and can handle cases where we're removing the head node*/
void list_remove(List* l, ListNode* n){

	/* if this node is both the front and back of the list, then we're going to end up with the empty list*/
	if (n == l->tail && n == l->tail){
		l->head = EMPTY_NODE;
		l->tail = EMPTY_NODE;
	} 
	/*else if this is the head of the list, but not the tail, we adjust the head*/
	else if (n == l->head){
		l->head = n->next;
	}
	/*else if this is the tail of the list, but not the head, we adjust the tail*/
	else if (n == l->tail){
		l->tail = n->last;
	}
	/*else this is neither the head, nor the tail, so we stitch up its two neighbour*/
	else {
		n->last->next = n->next;
		n->next->last = n->last;
	}
	
	/*now we reclaim some memory*/
	free(n);

}

void list_concat(List* dest, List* source){

	// if the source is empty
	if (source->head == EMPTY_NODE) return;


	if (dest->tail == EMPTY_NODE){
		// if the destination is empty, simply take source's list
		dest->head = source->head;
		dest->tail = source->tail;
	} else {
		// if the destination isn't empty, merge the two lists together
		dest->tail->next = source->head;	
		source->head->last = dest->tail;
		dest->tail = source->tail;
		
	}
	
	dest->size += source->size;

}

void list_delete(List* l){

	free(l);

}

void list_delete_nodes(List* l){

	ListNode* node;
	ListNode* next;

	node = l->head;
	
	/*keep going until we reach the end of the list*/
	while(true) {
		
		/*get the next node before we free this one*/
		next = node->next;
		
		/*free the node*/
		free(node);
		
		/*break if the next node is the empty list, otherwise assign next to node*/
		if (next == EMPTY_NODE) break;
		node = next;
	
	}
	
	list_delete(l);

}
