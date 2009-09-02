/*
  reorderRobust.c:
    Given a root node, reorder a tree either as postorder or preorder.
  Works on any valid tree, including those with singleton nodes and/or
  polytomies. Function inputs are derived from a phylo4 edge matrix. The
  new descendant node ordering is stored in descendantNew.
*/

#include <R.h>

typedef struct {
    int *descendantNew;
    int *ancestor;
    int *descendant;
    int nEdges;
    int index;
    } tree;

void postorderRobust(tree*, int node);
void preorderRobust(tree*, int node);

void reorderRobust(int *descendantNew, int *root, int *ancestor,
    int *descendant, int *nEdges, int *order) {

    tree tr;
    tr.ancestor = ancestor;
    tr.descendant = descendant;
    tr.descendantNew = descendantNew;
    tr.nEdges = *nEdges;
    tr.index = 0;

    if (*order==0) {
        postorderRobust(&tr, *root);
    } else if (*order==1) {
        preorderRobust(&tr, *root);
    } else {
        error("invalid order type");
    }

}

// postorder: continue traversing to the end, then record node
void postorderRobust(tree *tr, int node) {
    for (int i=0; i<tr->nEdges; i++) {
        if (tr->ancestor[i]==node) {
            postorderRobust(tr, tr->descendant[i]);
        }
    }
    tr->descendantNew[tr->index] = node;
    tr->index += 1;
}

// preorder: record node before continuing traversal
void preorderRobust(tree *tr, int node) {
    tr->descendantNew[tr->index] = node;
    tr->index += 1;
    for (int i=0; i<tr->nEdges; i++) {
        if (tr->ancestor[i]==node) {
            preorderRobust(tr, tr->descendant[i]);
        }
    }
}
