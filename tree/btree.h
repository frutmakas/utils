/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/01/13 01:25:29 $
 * $Revision: 1.1.2.1 $
 * $Id: btree.h,v 1.1.2.1 2004/01/13 01:25:29 syed Exp $
 ********************************************************************/

#ifndef _BINARY_TreE_H_
#define _BINARY_TreE_H_
#ifdef DOOOONNNNEEEEUUUUHHH
template<class T> class btree {
	public :
		btree *parent, *left, *right;
		T obj;
		btree(btree *Parent=NULL) {
			this->parent = Parent;
			left = right = 0;
		}
		btree(T Obj) {
			parent = left = right =0;
			this->obj = obj;
		}

		~btree(){
		}

		void InsertLeft(btree *Parent) {
			if(!Parent) throw CBTreeException(PARENT_NULL);
			if(!parent) throw CBTreeException(CHILD_PARENT_NOT_NULL);
			if(Parent->left==NULL) {
				Parent->left = this;
				this->parent = Parent;
			} else {
				this->left = Parent->left;
				Parent->left->parent = this;
				parent->left=this;
			}
		}

		void InsertRight(btree *Parent) {
			if(!Parent) throw CBTreeException(PARENT_NULL);
			if(!parent) throw CBTreeException(CHILD_PARENT_NOT_NULL);
			if(Parent->right==NULL) {
				Parent->right = this;
				this->parent = Parent;
			} else {
				this->right = Parent->right;
				Parent->right->parent = this;
				parent->right=this;

			}
		}

		void DestroyMeAndMyFamily(btree *me) {
			if(me->left) {
				btree::DestroyMeAndFamily(me->left);
				me->left = NULL;
			} 
			if(me->right) {
				btree::DestroyMeAndFamily(me->right);
				me->right = NULL;
			}
			delete me;
		}

		int GetAllObjInLevelX(int X, const TStack<T> &results) {
			if(X!=0 && left) {
				GetAllObjInLevelX(X-1, results)
			}
			// on arrive a X = 0 ou plus de left
			if(X==0) result
		
		}


};
#endif // donnnneeeeeeeuuuhhhh

#endif // ifndef