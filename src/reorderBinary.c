/*
  reorderBinary.c:
    Given a root node, reorder a tree either as postorder or preorder.
  Works only on binary trees, in which each internal node has exactly 2
  descendants. Function inputs are derived from a phylo4 edge matrix.
  The new descendant node ordering is stored in descendantNew.
*/

#include <R.h>

typedef struct {
    int *descendantNew;
    int *ancestor;
    int *left;
    int *right;
    int nEdges;
    int index;
    } tree;

void postorderBinary(tree*, int node);
void preorderBinary(tree*, int node);

void reorderBinary(int *descendantNew, int *root, int *ancestor, int *left,
    int *right, int *nEdges, int *order) {

    tree tr;
    tr.ancestor = ancestor;
    tr.left = left;
    tr.right = right;
    tr.descendantNew = descendantNew;
    tr.nEdges = *nEdges;
    tr.index = 0;

    if (*order==0) {
      postorderBinary(&tr, *root);
    } else if (*order==1) {
      preorderBinary(&tr, *root);
    } else {
      error("invalid order type");
    }

}

// postorder: continue traversing to the end, then record node
void postorderBinary(tree *tr, int node) {
    for (int i=0; i<tr->nEdges; i++) {
        if (tr->ancestor[i]==node) {
            postorderBinary(tr, tr->left[i]);
            postorderBinary(tr, tr->right[i]);
        }
     }
    tr->descendantNew[tr->index] = node;
    tr->index += 1;
}

// preorder: record node first, then continue traversing
void preorderBinary(tree *tr, int node) {
    tr->descendantNew[tr->index] = node;
    tr->index += 1;
    for (int i=0; i<tr->nEdges; i++) {
        if (tr->ancestor[i]==node) {
            preorderBinary(tr, tr->left[i]);
            preorderBinary(tr, tr->right[i]);
        }
    }
}
