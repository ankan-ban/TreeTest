// experiments with alpha-beta search on random tree

#include <stdio.h>
#include <stdlib.h>    
#include <math.h>
#include <time.h>
#include <assert.h>

// for timing CPU code : start
#include <windows.h>
double gTime;
LARGE_INTEGER freq;
#define START_TIMER { \
    LARGE_INTEGER count1, count2; \
    QueryPerformanceFrequency (&freq);  \
    QueryPerformanceCounter(&count1);

#define STOP_TIMER \
    QueryPerformanceCounter(&count2); \
    gTime = ((double)(count2.QuadPart-count1.QuadPart)*1000.0)/freq.QuadPart; \
    }
// for timing CPU code : end

//#define MAX_CHILDREN 40
#define MAX_CHILDREN 12
#define INF 10000.0f
const int g_depth = 10;

#define MAX_DEPTH 20


#define PV_NODE  1
#define CUT_NODE 2
#define ALL_NODE 3

struct Node
{
    float nodeVal;      // value from eval function for leaves, best searched value for interior nodes
    Node *children;     // pointer to array containing all child nodes
    Node *parent;       // pointer to parent node
    Node *best;         // pointer to best child

    unsigned char nChildren;       // no of child nodes
    unsigned char bestChild;       // most promising child/branch from this node
    unsigned char nodeType;        // PV, CUT or ALL node
    unsigned char nChildsExplored; // num of chlidren explored (only valid for CUT nodes)

    bool          isMaxNode;            // totally redundant, kept here for simplicity.
    int           frontierOffset;       // offset of first leaf in the frontier of the subtree whose root is this node
    int           numChildrenAtFrontier;// no of children of the subtree at frontier
};

int gTotalNodes;
int gLeafNodes;

void freeTree(Node *root)
{
    if (!root->children)
        return;

    for (int i = 0; i < root->nChildren;i++)
        freeTree(&root->children[i]);

    free(root->children);
}

void genTree(Node *root, int depth)
{
    gTotalNodes++;

    root->frontierOffset = -1;
    root->numChildrenAtFrontier = 0;


    if (g_depth % 2 == 0)
    {
        // even depths are maximizing nodes
        root->isMaxNode = !(depth % 2);
    }
    else
    {
        root->isMaxNode = !!(depth % 2);
    }

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
        children[i].nChildsExplored = 0;
        children[i].nChildren = 0;
        children[i].children = NULL;
        children[i].nodeVal = 0.0f;
        children[i].parent = root;

        genTree (&children[i], depth - 1);
    }

    root->nChildren = nChildren;
    root->children = children;
    root->nChildsExplored = 0;
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

bool isBetter(Node *node, float val)
{
    if (node->isMaxNode && val > node->nodeVal)
        return true;

    if ((!node->isMaxNode) && (val < node->nodeVal))
        return true;

    return false;
}

// returns true if val1 is better than val2
bool isBetter(bool isMaxNode, float val1, float val2)
{
    if (isMaxNode && (val1 > val2))
        return true;

    if ((!isMaxNode) && (val1 < val2))
        return true;

    return false;
}

int propogateFrontierOffsets(Node *node, int *childrenAtFrontier)
{
    if (node->frontierOffset != -1)
    {
        *childrenAtFrontier = node->numChildrenAtFrontier;
        return node->frontierOffset;
    }

    assert(node->nChildren > 0);
    int last, count;
    int first = propogateFrontierOffsets(&node->children[0], &count);
    if (node->nChildsExplored > 1)
    {
        for (int i = 1; i < node->nChildsExplored; i++)
        {
            last = propogateFrontierOffsets(&node->children[i], &count);
        }
        last += count;
        count = last - first;
    }

    node->frontierOffset = first;
    node->numChildrenAtFrontier = count;

    *childrenAtFrontier = count;

    return first;
}

int exploreSubTreeCount = 0;

bool expandNode(Node **fullCurrentFrontier, float *currentNodeVals, int i, float curBest, Node *subTreeRoot, bool *ignored);

// starts exploring a node as if it's a CUT node with cutVal as the value to check against
// returns either cutVal if a  better value couldn't be found,  or value of the best found node otherwise
// works only on CUT and ALL nodes, no node is marked as PV node by this function
float exploreSubTree(Node *node, float cutVal)
{
    exploreSubTreeCount++;

    // isMaxLevel is true if node's *parent* is a MAX level
    bool isMaxLevel = !node->isMaxNode;
    if (node->children == NULL)
    {
        // leaf node
        return isMaxLevel ? max(cutVal, node->nodeVal) 
                          : min(cutVal, node->nodeVal);
    }

    // full frontier
    Node **fullCurrentFrontier = NULL;
    Node **fullNextFrontier = NULL;

    // only for the last level
    float *currentNodeVals;

    int nCurr, nNext;

    // treat passed-in node as CUT node
    node->nodeType = CUT_NODE;
    fullCurrentFrontier = &node;
    nCurr = 1;

    bool currentIsMaxLevel = !isMaxLevel;
    int subDepth = 0;
    while (1)
    {
        bool secondLastLevel = false;
        secondLastLevel = (fullCurrentFrontier[0]->children[0].children == NULL);

        // explore the frontier nodes

        // figure out total no. of childs that need to be explored for all frontier nodes
        nNext = 0;
        for (int j=0; j<nCurr; j++)
        {
            switch (fullCurrentFrontier[j]->nodeType)
            {
                case ALL_NODE:
                    // if it's the second last level, we will explore all ALL nodes and find the best immediately
                    nNext += secondLastLevel ? 1 : fullCurrentFrontier[j]->nChildren;
                    break;
                case CUT_NODE:
                    nNext++;
                    break;
            }
        }
        
        // no more child nodes, we are at leaves..
        if (nNext == 0)
            break;

        // allocate memory for the next frontier
        fullNextFrontier = (Node**) malloc (nNext * sizeof(Node *));

        if (secondLastLevel) {
            assert(nCurr == nNext);
            currentNodeVals = (float *)malloc(nNext * sizeof(float));
        }
        
        int index = 0;
        // fill in the next frontier
        for (int j=0; j<nCurr; j++)
        {
            switch (fullCurrentFrontier[j]->nodeType)
            {
                case ALL_NODE:
                    if (secondLastLevel)
                    {
                        Node * curNode = fullCurrentFrontier[j];
                        curNode->nodeVal = curNode->isMaxNode ? -INF : INF;
                        for (int k = 0; k < curNode->nChildren; k++)
                        {
                            if (curNode->isMaxNode)
                            {
                                if (curNode->children[k].nodeVal > curNode->nodeVal)
                                {
                                    curNode->nodeVal = curNode->children[k].nodeVal;
                                    curNode->best = &curNode->children[k];
                                    curNode->bestChild = k;
                                }
                            }
                            else
                            {
                                if (curNode->children[k].nodeVal < curNode->nodeVal)
                                {
                                    curNode->nodeVal = curNode->children[k].nodeVal;
                                    curNode->best = &curNode->children[k];
                                    curNode->bestChild = k;
                                }
                            }
                        }
                        if (secondLastLevel)
                            currentNodeVals[index] = curNode->nodeVal;
                        fullNextFrontier[index++] = curNode;
                        curNode->nChildsExplored = curNode->nChildren;
                    }
                    else
                    {
                        for (int k = 0; k < fullCurrentFrontier[j]->nChildren; k++)
                        {
                            Node *curNode = &(fullCurrentFrontier[j]->children[k]);
                            curNode->nodeType = CUT_NODE;   // child of ALL node is cut node
                            fullNextFrontier[index++] = curNode;
                        }
                        fullCurrentFrontier[j]->nChildsExplored = fullCurrentFrontier[j]->nChildren;
                        fullCurrentFrontier[j]->best = &(fullCurrentFrontier[j]->children[0]);
                    }
                    break;
                case CUT_NODE:
                    {
                        fullCurrentFrontier[j]->nChildsExplored = 1;
                        fullCurrentFrontier[j]->best = &(fullCurrentFrontier[j]->children[0]);
                        fullCurrentFrontier[j]->children[0].nodeType = ALL_NODE;
                        fullNextFrontier[index] = &(fullCurrentFrontier[j]->children[0]);
                        if (secondLastLevel)
                        {
                            currentNodeVals[index] = fullNextFrontier[index]->nodeVal;
                            fullCurrentFrontier[j]->nodeVal = currentNodeVals[index];
                        }
                        index++;
                    }
            }

            if (secondLastLevel)
            {
                fullNextFrontier[j]->numChildrenAtFrontier = 1;
                fullNextFrontier[j]->frontierOffset = j;
            }
        }

        // go to next depth, set current = next
        if (subDepth != 0)
            free (fullCurrentFrontier);

        nCurr = nNext;
        fullCurrentFrontier = fullNextFrontier;

        subDepth++;

        currentIsMaxLevel = !currentIsMaxLevel;

        if (secondLastLevel)
            break;
    }

    int count;
    int start = propogateFrontierOffsets(node, &count);
    assert(start == 0 && count == nCurr);
    
    float *minScan = (float *)malloc(sizeof(float) * nCurr);
    float *maxScan = (float *)malloc(sizeof(float) * nCurr);
    bool *ignored = (bool *)malloc(sizeof(bool) * nCurr);
    memset(ignored, 0, sizeof(bool) * nCurr);

    float curMin, curMax;
    curMin = curMax = cutVal;  // init. with cut val

    float computedVal = isMaxLevel ? -INF : INF;

    int nRejected;
    int nExpnded;

    do
    {
        nRejected = 0;
        nExpnded = 0;

        curMin = curMax = cutVal;

        for (int i = 0; i < nCurr; i++)
        {
            // all nodes that are explored here must be ALL nodes
            assert(fullCurrentFrontier[i]->nodeType = ALL_NODE ||
                fullCurrentFrontier[i]->nChildsExplored == fullCurrentFrontier[i]->nChildren);

            if (ignored[i])
                continue;

            minScan[i] = curMin;
            maxScan[i] = curMax;

            // this node was expected to be less than the PV node value
            if (isMaxLevel)
            {
                if (currentNodeVals[i] <= curMin)
                {
                    // reject this branch (i.e, no need to evaluate any more siblings)
                    nRejected++;
                    ignored[i] = true;
                }
                if (currentNodeVals[i] > curMax)
                {
                    curMax = currentNodeVals[i];
                    // need to evaluate more siblings of this node
                    if (expandNode(fullCurrentFrontier, currentNodeVals, i, curMax, node, ignored))
                        nExpnded++;
                    computedVal = currentNodeVals[i];
                }
            }
            else //if (!isMaxLevel)
            {
                if (currentNodeVals[i] >= curMax)
                {
                    // reject this branch (i.e, no need to evaluate any more siblings of this node)
                    nRejected++;
                    ignored[i] = true;
                }
                if (currentNodeVals[i] < curMin)
                {
                    curMin = currentNodeVals[i];
                    // need to evaluate more siblings of this node
                    if(expandNode(fullCurrentFrontier, currentNodeVals, i, curMin, node, ignored))
                        nExpnded++;
                    computedVal = currentNodeVals[i];
                }
            }
        }

        // the loop is done when nothing gets expanded anymore
    } while (nExpnded);

    return isMaxLevel ? max(cutVal, computedVal)
                      : min(cutVal, computedVal);
}

// returns true if the node was actually expanded (sibling evaluated)
//         false otherwise (if there are no siblings, or if the node is a PV or subTreeRoot)
bool expandNode(Node **fullCurrentFrontier, float *currentNodeVals, int i, float curBest, Node *subTreeRoot, bool *ignored)
{
    Node *thisNode = fullCurrentFrontier[i];
    bool isMaxLevel = thisNode->isMaxNode;

    // no need to explore any siblings if the current node itself is subtree root!
    if (thisNode == subTreeRoot)
        return false;

    Node *parentNode = thisNode->parent;

    assert(parentNode->nodeType == CUT_NODE);

    Node *currentParent = parentNode;
    while (currentParent->nChildsExplored == currentParent->nChildren)
    {
        // break out early if we reached the subTreeRoot!
        if (currentParent == subTreeRoot)
        {
            //if (isBetter(currentParent, currentNodeVals[i]))
            {
                currentParent->best = thisNode;
                currentParent->nodeVal = currentNodeVals[i];
                thisNode->nodeType = PV_NODE;
            }
            return false;
        }
        // already explored the last child, need to explore more nodes of GrandParent of Parent node?
        // parent of parent should be an ALL node, with everything explored already ?
        currentParent->nodeVal = currentNodeVals[i];
        thisNode = currentParent;
        currentParent = currentParent->parent;

        if (currentParent->nodeType == ALL_NODE)
        {
            //if (isBetter(currentParent, currentNodeVals[i]))

            // when moving up, only cross an ALL node only when 
            //  i)  either it's the last child
            //  ii) there is no 'open' child towards the right of current node
            //     - open CHILD are the nodes that have possibility of being better than the current node
            bool isBest = true;
            // bool onRight = false;

            // TODO: Check only nodes towards the right. No need to check all nodes of parent
            // ANKAN - BUG! TODO!!! just checking nodes towards right doesn't work. Figure out why?

            for (int n = 0; n < currentParent->nChildren; n++)
            {
                if ((&currentParent->children[n]) != thisNode)
                {
                    //if (!onRight)
                    //    continue;

                    for (int k = 0; k < currentParent->children[n].numChildrenAtFrontier; k++)
                    {
                        int frontierOffset = currentParent->children[n].frontierOffset + k;
                        if ((!ignored[frontierOffset]) &&
                            isBetter(currentParent->isMaxNode, currentNodeVals[frontierOffset], currentNodeVals[i]))
                        {
                            isBest = false;
                        }
                        else
                        {
                            ignored[frontierOffset] = true;
                        }
                    }
                    if (!isBest)
                    {
                        break;
                    }
                }
                //else
                //{
                //    onRight = true;
                //}
            }


            if (isBest)
            {
                currentParent->best = thisNode;
                currentParent->nodeVal = currentNodeVals[i];
            }
            else
            {
                // wait for other better nodes to be sorted, out
                // don't move up
                return false;
            }
        }

        if (currentParent->nodeType == PV_NODE)
        {
            // Ankan TODO: might need to propogate up only if it's actually better
            currentParent->best = thisNode;
            thisNode->nodeType = PV_NODE;


            // currentParent->nodeVal = currentNodeVals[i];

            // need to propogate the value all the way to root!
            while (currentParent)
            {
                currentParent->nodeVal = currentNodeVals[i];

                // stop when the chain of PV nodes breaks!
                if (currentParent->nodeType != PV_NODE)
                    break;
                currentParent = currentParent->parent;
            }

            return false;
        }
    }

    assert(currentParent->nodeType != PV_NODE);
    {
        // here we have a parent node with unexplored children, explore one of them
        //assert(currentParent->nodeType == CUT_NODE);
        Node *sibling = &(currentParent->children[currentParent->nChildsExplored]);
        sibling->nodeType = ALL_NODE;
        currentParent->nChildsExplored++;
        float siblingVal = exploreSubTree(sibling, curBest);
        if (currentParent->isMaxNode && siblingVal > currentNodeVals[i])
        {
            currentNodeVals[i] = siblingVal;
            currentParent->best = sibling;
            fullCurrentFrontier[i] = sibling;
        }
        if ((!currentParent->isMaxNode) && siblingVal < currentNodeVals[i])
        {
            currentNodeVals[i] = siblingVal;
            currentParent->best = sibling;
            fullCurrentFrontier[i] = sibling;
        }

    }
    currentParent->nodeVal = currentNodeVals[i];

    // propogate the value all the way up if the parent's value was coming from this node
    Node *cur = currentParent;
    while (true)
    {
        Node *par = cur->parent;
        if (par->best == cur)
        {
            par->nodeVal = cur->nodeVal;
            cur = par;
        }
        else
        {
            break;
        }
    }

    return true;
}

// a non-recursive and hopefully somewhat parallel algorithm based on alpha beta
float exploreTree(Node *node, int depth)
{
    exploreSubTreeCount = 0;

    // the frontier / current list of nodes that need to be explored / nodes at the current level
    // TODO: many of these lists are probably redundant - get rid of some later
    Node *currentPVNode   = NULL;       Node *nextPVNode   = NULL;
    Node *currentCutNodes = NULL;       Node *nextCutNodes = NULL;
    Node *currentAllNodes = NULL;       Node *nextAllNodes = NULL;
    int nCutNodes = 0, nAllNodes = 0;
    int nNextCutNodes = 0, nNextAllNodes = 0;

    // full frontier
    Node **fullCurrentFrontier = NULL;
    Node **fullNextFrontier    = NULL;
    int nCurr, nNext;

    // 1. PV based initial tree generation
    node->nodeType = PV_NODE;   // the root is always a PV node
    currentPVNode = node;
    fullCurrentFrontier = &node;
    nCurr = 1;

    bool isMaxLevel = true;
    for (int i = 0; i < depth - 1; i++)
    {
        // explore the frontier nodes

        // figure out total no. of childs that need to be explored for all frontier nodes
        nNext = 0;
        for (int j=0; j<nCurr; j++)
        {
            int childsToExplore = 0;
            switch (fullCurrentFrontier[j]->nodeType)
            {
                case PV_NODE:
                case ALL_NODE:
                    childsToExplore = fullCurrentFrontier[j]->nChildren;
                    break;
                case CUT_NODE:
                    childsToExplore = 1;
                    break;
            }

            nNext += childsToExplore;
        }
        
        // allocate memory for the next frontier
        fullNextFrontier = (Node**) malloc (nNext * sizeof(Node *));
        
        int index = 0;
        // fill in the next frontier
        for (int j=0; j<nCurr; j++)
        {
            switch (fullCurrentFrontier[j]->nodeType)
            {
                case PV_NODE:
                    for (int k = 0; k < fullCurrentFrontier[j]->nChildren; k++)
                    {
                        Node *curNode = &(fullCurrentFrontier[j]->children[k]);
                        if (index == 0) // first child of PV node is PV node
                            curNode->nodeType = PV_NODE;
                        else            // others are CUT nodes
                            curNode->nodeType = CUT_NODE;
                        fullNextFrontier[index++] = curNode;
                    }
                    fullCurrentFrontier[j]->nChildsExplored = fullCurrentFrontier[j]->nChildren;

                    break;
                case ALL_NODE:
                    for (int k = 0; k < fullCurrentFrontier[j]->nChildren; k++)
                    {
                        Node *curNode = &(fullCurrentFrontier[j]->children[k]);
                        curNode->nodeType = CUT_NODE;   // child of ALL node is cut node
                        fullNextFrontier[index++] = curNode;
                    }
                    fullCurrentFrontier[j]->nChildsExplored = fullCurrentFrontier[j]->nChildren;
                    break;
                case CUT_NODE:
                    {
                        fullCurrentFrontier[j]->nChildsExplored = 1;
                        fullCurrentFrontier[j]->children[0].nodeType = ALL_NODE;
                        fullNextFrontier[index++] = &(fullCurrentFrontier[j]->children[0]);
                    }
            }

            // Ankan - not known yet, but initialize with first child
            fullCurrentFrontier[j]->best = &(fullCurrentFrontier[j]->children[0]);
            fullCurrentFrontier[j]->bestChild = 0;
        }

        assert(index == nNext);

        // go to next depth, set current = next
        if (i!=0)
            free (fullCurrentFrontier);

        nCurr = nNext;
        fullCurrentFrontier = fullNextFrontier;

        isMaxLevel = !isMaxLevel;
    }
    

    // when generating the last level evaluate all ALL nodes at the level just above the CUT node leaves
    // for depth 5 search, we need to do a MAX reduction (see modern gpu's segmented reduction example when implementing parallel version)

    fullNextFrontier = (Node**) malloc (nCurr * sizeof(Node *));
    bool *expectedMore = (bool*) malloc (nCurr * sizeof(bool));
    float *currentNodeVals = (float *) malloc(nCurr * sizeof(float));

    for (int i=0;i<nCurr;i++)
    {
        Node *curNode = fullCurrentFrontier[i];
        curNode->nodeVal = isMaxLevel ? -INF : INF;
        if(curNode->nodeType == PV_NODE || curNode->nodeType == ALL_NODE)
        {
            for (int k = 0; k < curNode->nChildren; k++)
            {
                if (isMaxLevel)
                {
                    if (curNode->children[k].nodeVal > curNode->nodeVal)
                    {
                        curNode->nodeVal = curNode->children[k].nodeVal;
                        curNode->best = &curNode->children[k];
                        curNode->bestChild = k;
                    }
                }
                else
                {
                    if (curNode->children[k].nodeVal < curNode->nodeVal)
                    {
                        curNode->nodeVal = curNode->children[k].nodeVal;
                        curNode->best = &curNode->children[k];
                        curNode->bestChild = k;
                    }
                }
            }
            fullNextFrontier[i] = curNode;
            currentNodeVals[i] = curNode->nodeVal;
            curNode->nChildsExplored = curNode->nChildren;
            expectedMore[i] = isMaxLevel ? false : true;
        }
        else
        {
            curNode->best = &curNode->children[0];
            curNode->bestChild = 0;

            curNode->nChildsExplored = 1;
            curNode->children[0].nodeType = ALL_NODE;
            fullNextFrontier[i] = &curNode->children[0];
            currentNodeVals[i] = fullNextFrontier[i]->nodeVal;
            expectedMore[i] = isMaxLevel ? true : false;
        }

        fullNextFrontier[i]->numChildrenAtFrontier = 1;
        fullNextFrontier[i]->frontierOffset = i;
    }

    free (fullCurrentFrontier);
    fullCurrentFrontier = fullNextFrontier;


    float *minScan = (float *) malloc (sizeof(float) * nCurr);
    float *maxScan = (float *) malloc (sizeof(float) * nCurr);
    bool *ignored = (bool *)malloc(sizeof(bool) * nCurr);
    memset(ignored, 0, sizeof(bool) * nCurr);

    float curMin, curMax;
    curMin = curMax = currentNodeVals[0];  // init. with value of PV node

    int nRejected = 0;
    int nExpnded = 0;

    // propogate frontier offsets from leafs up the tree
    int count;
    int start = propogateFrontierOffsets(node, &count);
    assert(start == 0 && count == nCurr);

    int iterations = 0;

    do
    {
        iterations++;

        nRejected = 0;
        nExpnded = 0;

        // ignore all nodes to the left of bestNow
        for (int i = 0; i < nCurr; i++)
        {
            ignored[i] = true;

            // Unnecessary and expensive. TODO: figure out better way
            bool found = false;
            Node *bestNow = node;
            while (bestNow->children)
            {
                bestNow = bestNow->best;
                if (bestNow == fullCurrentFrontier[i])
                {
                    found = true;
                }
            }

            if (found)
            {
                curMin = curMax = bestNow->nodeVal;
                break;
            }
        }

        for (int i = 1; i < nCurr; i++)
        {
            // all nodes that are explored here must be ALL nodes
            assert(fullCurrentFrontier[i]->nodeType = ALL_NODE ||
                   fullCurrentFrontier[i]->nChildsExplored == fullCurrentFrontier[i]->nChildren);

            if (ignored[i])
                continue;

            minScan[i] = curMin;
            maxScan[i] = curMax;

            // this node was expected to be less than the PV node value
            if (expectedMore[i] == false)
            {
                if (currentNodeVals[i] <= curMin)
                {
                    // reject this branch (i.e, no need to evaluate any more siblings)
                    nRejected++;
                    ignored[i] = true;
                }
                if (currentNodeVals[i] > curMax)
                {
                    curMax = currentNodeVals[i];
                    // need to evaluate more siblings of this node
                    expandNode(fullCurrentFrontier, currentNodeVals, i, curMax, NULL, ignored);
                    nExpnded++;
                }
            }
            else //if (expectedMore[i])
            {
                if (currentNodeVals[i] >= curMax)
                {
                    // reject this branch (i.e, no need to evaluate any more siblings of this node)
                    nRejected++;
                    ignored[i] = true;
                }
                if (currentNodeVals[i] < curMin)
                {
                    curMin = currentNodeVals[i];
                    // need to evaluate more siblings of this node
                    expandNode(fullCurrentFrontier, currentNodeVals, i, curMin, NULL, ignored);
                    nExpnded++;
                }
            }
        }

    } while (nExpnded);

    printf("\nFrontier Nodes: %d, main loop iterations: %d, explore subtree count: %d", nCurr, iterations, exploreSubTreeCount);
    printf("\nExplore Tree found value: %f\n", node->nodeVal);



    free (fullCurrentFrontier);
    return node->nodeVal;
}


struct ListItem
{
    Node *node;
    float merit;    // upper bound for live, actual value for solved
    int   depth;    // depth of the node (used to identify leaf/root)
    bool  live;     // true - live, false - solved
};

#define MAX_LIST_SIZE 100000


int g_sssNodes = 0;
class List
{
private:
    // simple array based implementation. O(n) complexity for most operations
    ListItem m_list[MAX_LIST_SIZE];
    int n;

public:
    List()
    {
        n = 0;
        memset(m_list, 0, sizeof(m_list));
    };

    void addItem(Node *node, bool live, float merit, int depth)
    {
        ListItem newItem;
        newItem.live = live;
        newItem.node = node;
        newItem.merit = merit;
        newItem.depth = depth;
        m_list[n++] = newItem;

        if (live == true)
            g_sssNodes++;
    };


    // return index of max node in the list
    int findMax()
    {
        int maxIndex = 0;
        ListItem max = m_list[0];
        for (int i=1; i<n; i++)
        {
            if (m_list[i].merit > max.merit)
            {
                max = m_list[i];
                maxIndex = i;
            }
        }
        
        return maxIndex;
    }

    ListItem getMax()
    {
        return m_list[findMax()];
    };

    // return index of node in the list
    int findItem(Node *node)
    {
        for (int i=0; i<n; i++)
            if (m_list[i].node == node)
                return i;

        // NOT FOUND
        return -1;
    }
    
    void deleteIndex(int index)
    {
        if (index >= 0)
        {
            for (int i=index; i < n-1; i++)
                m_list[i] = m_list[i+1];
            n--;
        }
    }

    void deleteItem(Node *node)
    {
        deleteIndex(findItem(node));
    }

    ListItem extractMax()
    {
        int maxIndex = findMax();
        ListItem max = m_list[maxIndex];
        deleteIndex(maxIndex);
        return max;
    }

};

void purgeSubTree(List *list, Node *node)
{
    if (node->nChildsExplored)
    {
        for (int i=0; i<node->nChildsExplored; i++)
        {
            purgeSubTree(list, &node->children[i]);
        }
    }

    node->nChildsExplored = 0;

    list->deleteItem(node);
}

float SSS_star(Node *node, int depth)
{
    List *activeNodes = new List();

    node->nChildsExplored = 0;
    activeNodes->addItem(node, true, INF, 0);
    
    while(true)
    {
        ListItem node = activeNodes->extractMax();
        if (node.live)
        {
            if (node.depth == depth)    // leaf
            {
                activeNodes->addItem(node.node, false, min(node.merit, node.node->nodeVal), node.depth);                
            }
            else if (node.depth % 2 == 1)   // min node
            {
                activeNodes->addItem(&node.node->children[0], true, node.merit, node.depth + 1);

                node.node->nChildsExplored = 1;
            }
            else    // max node
            {
                for (int j=0; j<node.node->nChildren; j++)
                    activeNodes->addItem(&node.node->children[j], true, node.merit, node.depth + 1);

                node.node->nChildsExplored = node.node->nChildren;
            }
        }
        else    // solved
        {
            if (node.depth == 0)
            {
                return node.merit;
            }
            else if (node.depth % 2 == 1)   // min node
            {
                // purge parent and it's all children present in the list
                // TODO: optimize this
                purgeSubTree(activeNodes, node.node->parent);
                activeNodes->addItem(node.node->parent, false, node.merit, node.depth - 1);

                if (node.depth == 1)
                {
                    // find the child node index for best move
                    for (int  i=0; i<node.node->parent->nChildren; i++)
                    {
                        if (node.node == &node.node->parent->children[i])
                        {
                            node.node->parent->bestChild = i;
                            break;
                        }
                    }
                }
            }
            else    // max node
            {
                if (node.node->parent->nChildsExplored != node.node->parent->nChildren)
                {   // if node has unexplored brother, explore it
                    activeNodes->addItem(&node.node->parent->children[node.node->parent->nChildsExplored++], true, node.merit, node.depth);
                }
                else
                {
                    activeNodes->addItem(node.node->parent, false, node.merit, node.depth - 1);
                }
            }
        }
    }

    delete activeNodes;
    return 0;
}


int main2()
{
    printf("\n\nSize of node is %zd bytes\n\n", sizeof(Node));
    int randSeed = time(NULL);
    printf("Random Seed: %d\n", randSeed);
    srand(randSeed);
    // 3, and 4 are good. 6,7,8 are very good
    // 10 is excellent
    srand(10);
    // generate a random tree
    printf ("generating random tree of depth %d\n", g_depth);
    Node root = {0};
    START_TIMER
    genTree(&root, g_depth);
    STOP_TIMER
    printf("random tree generated, total nodes: %d, leaf nodes: %d, time: %g ms\n", gTotalNodes, gLeafNodes, gTime);



    float bestVal = 0;
    // search the best move using min-max search
    printf("searching the tree using min-max\n");
    START_TIMER
    bestVal = negaMax(&root, g_depth, g_depth);
    STOP_TIMER
    printf ("best move %d, score: %f, time: %g\n", root.bestChild, root.nodeVal, gTime);

    // search the best move using alpha-beta search
    printf("searching the tree using alpha-beta\n");
    START_TIMER
    bestVal = alphabeta(&root, g_depth, g_depth, -INF, INF);
    STOP_TIMER
    printf ("best move %d, score: %f\nnodes visited (leaves/interior/total): %d/%d/%d\n", 
            root.bestChild, root.nodeVal, gLeafNodesVisited, gInteriorNodesVisited, gLeafNodesVisited + gInteriorNodesVisited);
    printf("time taken: %g\n", gTime);

    START_TIMER
    exploreTree(&root, g_depth);
    STOP_TIMER
    printf("time taken: %g\n", gTime);    
    

    float val;
    START_TIMER
    val = SSS_star(&root, g_depth);
    STOP_TIMER
    printf("SSS* best node: %d, score: %f, nodes explored: %d, time taken: %g\n", root.bestChild, val, g_sssNodes, gTime);

    freeTree(&root);
    getchar();

    return 0;
}


int main()
{
    for (int i = 0; i < 1000; i++)
    {
        int randSeed = i + time(NULL);

        gTotalNodes = 0;
        gLeafNodes = 0;
        gLeafNodesVisited = 0;
        gInteriorNodesVisited = 0;


        printf("Random Seed: %d\n", randSeed);
        srand(randSeed);

        printf("generating random tree of depth %d\n", g_depth);
        Node root = { 0 };
        START_TIMER
            genTree(&root, g_depth);
        STOP_TIMER
        printf("random tree generated, total nodes: %d, leaf nodes: %d, time: %g ms\n", gTotalNodes, gLeafNodes, gTime);

        float bestValAB = 0;
        float bestValET = 0;
        // search the best move using alpha-beta search
        printf("searching the tree using alpha-beta\n");
        START_TIMER
            bestValAB = alphabeta(&root, g_depth, g_depth, -INF, INF);
        STOP_TIMER
        printf("best move %d, score: %f\nnodes visited (leaves/interior/total): %d/%d/%d\n",
               root.bestChild, root.nodeVal, gLeafNodesVisited, gInteriorNodesVisited, gLeafNodesVisited + gInteriorNodesVisited);
        printf("time taken: %g\n", gTime);


        START_TIMER
            bestValET = exploreTree(&root, g_depth);
        STOP_TIMER
        printf("time taken: %g\n\n\n\n", gTime);

        if (bestValET != bestValAB)
        {
            printf("\n*Mismatch found!*\n");
            getchar();
        }

        freeTree(&root);
        //getchar();
    }
    getchar();

    return 0;
}