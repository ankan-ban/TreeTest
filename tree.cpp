// experiments with alpha-beta search on random tree

#include <stdio.h>
#include <stdlib.h>    
#include <math.h>

#define MAX_CHILDREN 64
#define INF 10000

struct Node
{
    float nodeVal;      // value from eval function for leaves, best searched value for interior nodes
    int bestChild;      // most promising child/branch from this node

    int nChildren;      // no of child nodes
    Node *children;     // pointer to array containing all child nodes

    Node *parent;       // pointer to parent node
};


int gTotalNodes;
int gLeafNodes;

void genTree(Node *root, int depth)
{
    gTotalNodes++;

    if (depth == 0)
    {
        gLeafNodes++;
        // generate random node val (between 0 and 100) and return
        root->nodeVal = (rand() % 10000) / 100.0f;
        return;
    }

    // allocate memory for random number of children (have at least one children for now - to be changed later!)
    int nChildren = rand() % MAX_CHILDREN + 1;

    Node *children = (Node *) malloc (nChildren * sizeof(Node));

    for (int i=0; i<nChildren; i++)
    {
        children[i].nChildren = 0;
        children[i].children = NULL;
        children[i].nodeVal = 0.0f;
        children[i].parent = root;

        genTree (&children[i], depth - 1);
    }

    root->nChildren = nChildren;
    root->children = children;
}


float negaMax(Node *node, int depth, int origDepth)
{
    if (depth == 0)
    {
        // eval for even depths, -eval for odd depths
        if (origDepth % 2 == 0)
            return node->nodeVal;
        else
            return -node->nodeVal;
    }

    // choose the best child
    float bestScore = -INF;
    int bestChild = 0;

    for (int i = 0; i < node->nChildren; i++)
    {
        float curScore = -negaMax(&node->children[i], depth - 1, origDepth);
        if (curScore > bestScore)
        {
            bestScore = curScore;
            bestChild = i;
        }
    }

    node->nodeVal = bestScore;
    node->bestChild = bestChild;

    return bestScore;
}



int gInteriorNodesVisited = 0;
int gLeafNodesVisited = 0;

float alphabeta(Node *node, int depth, int origDepth, float alpha, float beta)
{
    if (depth == 0)
    {
        gLeafNodesVisited++;

        // eval for even depths, -eval for odd depths
        if (origDepth % 2 == 0)
            return node->nodeVal;
        else
            return -node->nodeVal;
    }

    gInteriorNodesVisited++;

    // choose the best child
    int bestChild = 0;

    for (int i = 0; i < node->nChildren; i++)
    {
        float curScore = -alphabeta(&node->children[i], depth - 1, origDepth, -beta, -alpha);
        if (curScore >= beta)
        {
            return beta;
        }

        if (curScore > alpha)
        {
            alpha = curScore;
            bestChild = i;
        }
    }

    node->nodeVal = alpha;
    node->bestChild = bestChild;

    return alpha;

}

int main()
{
    int depth = 5;


    // generate a random tree
    printf ("generating random tree of depth %d\n", depth);
    Node root = {0};
    genTree(&root, depth);
    printf("random tree generated, total nodes: %d, leaf nodes: %d\n", gTotalNodes, gLeafNodes);

    float bestVal = 0;
    // search the best move using min-max search
    printf("searching the tree using min-max\n");
    bestVal = negaMax(&root, depth, depth);
    printf ("best move %d, score: %f\n", root.bestChild, root.nodeVal);

    // search the best move using alpha-beta search
    printf("searching the tree using alpha-beta\n");
    bestVal = alphabeta(&root, depth, depth, -INF, INF);
    printf ("best move %d, score: %f\nnodes visited (leaves/interior/total): %d/%d/%d\n", 
            root.bestChild, root.nodeVal, gLeafNodesVisited, gInteriorNodesVisited, gLeafNodesVisited + gInteriorNodesVisited);

    return 0;
}