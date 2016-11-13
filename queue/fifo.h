/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/30 04:23:01 $
 * $Revision: 1.1.2.6 $
 * $Id: fifo.h,v 1.1.2.6 2006/07/30 04:23:01 syed Exp $
 ********************************************************************/
 
#ifndef _FIFO_H_
#define _FIFO_H_

namespace utilitis {

template<class T> class TStackObj;

template<class T> class TStackObj {
	public:
		T obj; TStackObj<T> *next;

		TStackObj(const T& newobj){
			obj = newobj;
			next=NULL;
		}
		~TStackObj(){};
};

template<class T> class TStack;

template<class T> class TStack {
	private :
		TStackObj<T> *head, *tail;
        int items;
        bool locked;
	public:
		int Stack(const T & newobj){
            if(locked) return 0;
            locked=true;
			if (head==NULL) {
				head = new TStackObj<T>(newobj);
				head->next=NULL;
				tail = head;
                items++;
                locked = false;
				return 1;
			}

			if (tail==NULL && head!=NULL) {
				// y a un blem
                locked = false;
				return 0;
			}

			if (tail!=NULL) {
				tail->next = new TStackObj<T>(newobj);
				tail = tail->next;
                items++;
                locked = false;
				return 1;
			}
		}
		int Unstack(T &obj){
            if(locked) return 0;
            locked = true;
			if (head==NULL) {
                locked =false;
                return 0;
            }
			obj = head->obj;
			TStackObj<T> *tmp = head;
			if(head==tail) { // on a un seul element ...  head = tail = null
                head = tail = NULL;
            } else {
                head = head->next;
            }
			delete tmp;
            items--;
            locked = false;
			return 1;
		}
		T Spy(){ // read head without unstack
			if (head==NULL) return T();
			return head->obj;
		}

		TStack(){head = tail = NULL; items=0; locked=false;}
		~TStack(){
            locked = true;
			while (head) {
				TStackObj<T> *obj = head;
				head = head->next;
				delete obj;
			}
            head=tail=NULL;
            items=0;
            locked=false;
		}
		bool isEmpty() {return (head==NULL);}
        int Count() { return items; }

        int Flush() {
            if(locked) return 0;
            locked =true;
			while (head) {
				TStackObj<T> *obj = head;
				head = head->next;
				delete obj;
			}
            items=0;
            head=tail=NULL;
            locked = false;
            return 1;
        }
};
}
#endif



