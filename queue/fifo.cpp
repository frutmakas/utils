/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/04/05 18:49:16 $
 * $Revision: 1.1.2.4 $
 * $Id: fifo.cpp,v 1.1.2.4 2004/04/05 18:49:16 syed Exp $
 ********************************************************************/
 
#ifdef APESALLETAKJADI
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#endif
#include "queue/fifo.h"

template<class T> TStackObj<T>::TStackObj(const T& newobj) {
	obj = newobj;
	next=NULL;
}

template<class T> TStackObj<T>::~TStackObj(){
	delete obj;
}

template<class T> void TStack<T>::Stack(const T & newobj){
	if (head==NULL) {
		head = new TStackObj<T>(newobj);
		head->next=NULL;
		tail = head;
		return;
	}

	if (tail==NULL && head!=NULL) {
		// y a un blem
	}

	if (tail!=NULL) {
		tail->next = new TStackObj<T>(newobj);
		tail = tail->next;
		return;
	}
}

template<class T> T TStack<T>::Unstack(){
	if (head==NULL) return T();
	T obj = head->obj;
	TStackObj<T> *tmp = head;
	head = head->next;
	delete tmp;
	return obj;
}

template<class T> TStack<T>::TStack(){
	head = tail = NULL;
}
		
template<class T> TStack<T>::~TStack(){
	while (head) {
		TStackObj<T> *obj = head; 
		head = head->next;
		delete obj;
	}
}

template<class T> bool TStack<T>::isEmpty(){
	return (head==NULL);
}

#endif
