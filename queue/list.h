/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/04/05 18:49:17 $
 * $Revision: 1.1.2.9 $
 * $Id: list.h,v 1.1.2.9 2004/04/05 18:49:17 syed Exp $
 ********************************************************************/
 
#ifndef _List_LIST_H_
#define _List_LIST_H_

#define BUGGY
#ifdef BUGGY
#include <stdio.h>
#endif
#define OUT_OF_RANGE				1010
#define SORT_INTEGRITY_VIOLATION	1011

class CListException{
public:
	int err_code;
	CListException(int err){
		err_code=err;
	}
};

template<class T> class TSortedListObj;

template<class T> class TSortedListObj {
	public:
		T obj; TSortedListObj<T> *next, *prev;
		
		TSortedListObj(const T& newobj){
			obj = newobj;
			prev=next=NULL;
		}
		~TSortedListObj(){};
};

template<class T> class TSortedList;

template<class T> class TSortedList {
	private : 
		TSortedListObj<T> *head;
	public:
		int count;
		void Insert(const T & newobj){
			TSortedListObj<T> *tmp = head;
			count++;
			if (tmp==NULL) {
				head = new TSortedListObj<T>(newobj);
				return;
			}

			while(tmp->next!= NULL) {
				if (tmp->obj>newobj) {
					TSortedListObj<T> *casseur= new TSortedListObj<T>(newobj);
					if(tmp->prev==NULL) { //insert at first element
						casseur->next=head;
						head->prev=casseur;
						head=casseur;
						return;
					}
					casseur->prev = tmp->prev;
					casseur->next = tmp;
					tmp->prev = casseur;
					casseur->prev->next=casseur;
					return;
				}
				tmp=tmp->next;
			}
			// we are at the end of list
			if(tmp->obj < newobj) {
				TSortedListObj<T> *casseur= new TSortedListObj<T>(newobj);
				casseur->prev= tmp;
				tmp->next=casseur;
				return;
			}

			if(tmp->obj > newobj) { 
				TSortedListObj<T> *casseur = new TSortedListObj<T>(newobj);
				if(tmp->prev==NULL) {// at first element
					casseur->next = head;
					head->prev=casseur;
					head=casseur;
					return;
				}

				casseur->next = tmp;
				casseur->prev = tmp->prev;
				tmp->prev=casseur;
				casseur->prev->next=casseur;
				return;
			}
		}

		T Get(int element_number) const { // element_number starts with 1
			TSortedListObj<T> *tmp=head;
			while(tmp!=NULL) {
				if(--element_number==0) {
					return tmp->obj;
				}
				tmp=tmp->next;
			}
			throw CListException(OUT_OF_RANGE);
		}

		T Delete(int element_number){ // element_number starts with 1
			TSortedListObj<T> *tmp=head;
			while(tmp!=NULL) {
				if(--element_number==0) {
					T xobj = tmp->obj;
					if (tmp->prev == NULL) {
						head=tmp->next;
						head->prev=NULL;
						delete tmp;
						count--;
						return xobj;
					} 

					tmp->prev->next = tmp->next;

					if (tmp->next==NULL) {
						delete tmp;
						count--;
						return xobj;
					}

					tmp->next->prev = tmp->prev;
					delete tmp;
					count--;
					return xobj;
				}
				tmp=tmp->next;
			}
			throw CListException(OUT_OF_RANGE);
		}

		void Update(int element_number, T replacement_obj){
			if(element_number==0) throw CListException(OUT_OF_RANGE);
			TSortedListObj<T> *tmp=head;
			while(tmp!=NULL) {
				if(--element_number==0) {
					if(tmp->prev!=NULL){
						if(tmp->prev->obj>replacement_obj) throw CListException(SORT_INTEGRITY_VIOLATION);
					}
					if(tmp->next!=NULL){
						if(tmp->next->obj<replacement_obj) throw CListException(SORT_INTEGRITY_VIOLATION);
					}
					tmp->obj=replacement_obj;
					return;
				}
				tmp=tmp->next;
			}
			throw CListException(OUT_OF_RANGE);
		}

		int isExist(T search_obj) {
			int pos=0;
			TSortedListObj<T> *tmp=head;
			while(tmp!=NULL && tmp->obj <=search_obj) {
				pos++;
				if (tmp->obj==search_obj) {
					return pos;
				}
				tmp=tmp->next;
			}
			return -1;
		}

		int DeleteByValue(T search_obj) {
			int pos=0;
			TSortedListObj<T> *tmp=head;
			while(tmp!=NULL && tmp->obj <=search_obj) {
				pos++;
				if (tmp->obj==search_obj) {
					// delete here
					if (tmp->prev == NULL) {
						head=tmp->next;
						if (head != NULL) head->prev=NULL;
						delete tmp;
						count--;
						return 1;
					} 

					tmp->prev->next = tmp->next;

					if (tmp->next==NULL) {
						delete tmp;
						count--;
						return 1;
					}

					tmp->next->prev = tmp->prev;
					delete tmp;
					count--;
					return 1;
				}
				tmp=tmp->next;
			}
			return 0;
		}



		TSortedList(){head = NULL; count=0;}
		~TSortedList(){
			while (head) {
				TSortedListObj<T> *obj = head; 
				head = head->next;
				delete obj;
			}
			count=0;
		}

		TSortedList &operator=(const TSortedList &t){
			int i=count;
			TSortedListObj<T> *obj=NULL;
			for(i=1;i<=count;i++) {
				obj=head->next;
				delete head;
				head=obj;
			}

			T xobj;
			for(i=1;i<=t.count;i++) {
					xobj = t.Get(i);
					this->Insert(xobj);
			}
			return *this;
		}

		void copy(const TSortedList &t){
			int i=count;
			TSortedListObj<T> *obj=NULL;
			for(i=1;i<=count;i++) {
				obj=head->next;
				delete head;
				head=obj;
			}
			T xobj;
			for(i=1;i<=t.count;i++) {
					xobj = t.Get(i);
					this->Insert(xobj);
			}
			return *this;
		}


		bool isEmpty() {return (head==NULL);}
};

#endif
