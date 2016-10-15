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

#define MAX_CHILDREN 40
#define INF 10000

#define PV_NODE  1
#define CUT_NODE 2
#define ALL_NODE 3

struct Node
{
    float nodeVal;      // value from eval function for leaves, best searched value for interior nodes
    Node *children;     // pointer to array containing all child nodes
    Node *parent;       // pointer to parent node

    unsigned char nChildren;       // no of child nodes
    unsigned char bestChild;       // most promising child/branch from this node
    unsigned char nodeType;        // PV, CUT or ALL node
    unsigned char nChildsExplored; // num of chlidren explored (only valid for CUT nodes)
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

// a non-recursive and hopefully somewhat parallel algorithm based on alpha beta
float exploreTree(Node *node, int depth)
{
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
            switch (fullCurrentFrontier[j]->nodeType)
            {
                case PV_NODE:
                case ALL_NODE:
                    nNext += fullCurrentFrontier[j]->nChildren;
                    break;
                case CUT_NODE:
                    nNext++;
                    break;
            }
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
        }

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

    for (int i=0;i<nCurr;i++)
    {
        Node *curNode = fullCurrentFrontier[i];
        curNode->nodeVal = isMaxLevel ? -INF : INF;
        if(curNode->nodeType == PV_NODE || curNode->nodeType == ALL_NODE)
        {
            for (int k = 0; k < curNode->nChildren; k++)
            {
                curNode->nodeVal = isMaxLevel ? max(curNode->nodeVal, curNode->children[k].nodeVal) 
                                              : min(curNode->nodeVal, curNode->children[k].nodeVal);
            }
            fullNextFrontier[i] = curNode;
            curNode->nChildsExplored = curNode->nChildren;
            expectedMore[i] = isMaxLevel ? false : true;
        }
        else
        {
            curNode->nChildsExplored = 1;
            curNode->children[0].nodeType = ALL_NODE;
            fullNextFrontier[i] = &curNode->children[0];
            expectedMore[i] = isMaxLevel ? true : false;
        }
    }

    free (fullCurrentFrontier);
    fullCurrentFrontier = fullNextFrontier;


    float *minScan = (float *) malloc (sizeof(float) * nCurr);
    float *maxScan = (float *) malloc (sizeof(float) * nCurr);
    bool *ignored = (bool *)malloc(sizeof(bool) * nCurr);
    memset(ignored, 0, sizeof(bool) * nCurr);

    float curMin, curMax;
    curMin = curMax = fullCurrentFrontier[0]->nodeVal;  // init. with value of PV node

    int nRejected = 0;
    int nExpnded = 0;

    do
    {
        for (int i = 1; i < nCurr; i++)
        {
            // all nodes that are explored here must be ALL nodes
            assert(fullCurrentFrontier[i]->nodeType == ALL_NODE);

            if (ignored[i])
                continue;

            minScan[i] = curMin;
            maxScan[i] = curMax;

            // this node was expected to be less than the PV node value
            if (expectedMore[i] == false)
            {
                if (fullCurrentFrontier[i]->nodeVal <= curMin)
                {
                    // reject this branch (i.e, no need to evaluate any more siblings)
                    nRejected++;
                    ignored[i] = true;
                }
                if (fullCurrentFrontier[i]->nodeVal > curMax)
                {
                    curMax = fullCurrentFrontier[i]->nodeVal;
                    // need to evaluate more siblings of this node
                    
                    nExpnded++;
                }
            }
            else //if (expectedMore[i])
            {
                if (fullCurrentFrontier[i]->nodeVal >= curMax)
                {
                    // reject this branch (i.e, no need to evaluate any more siblings of this node)
                    nRejected++;
                    ignored[i] = true;
                }
                if (fullCurrentFrontier[i]->nodeVal < curMin)
                {
                    curMin = fullCurrentFrontier[i]->nodeVal;
                    // need to evaluate more siblings of this node
                    nExpnded++;
                }
            }
        }

    } while (nExpnded);

    // how to expand siblings of a node depends on the depth of the node
    // if it's a leaf, we can just generate the next sibling and be done with it.
    // ... at least unless is the last sibling.. if it's the last sibling.. need to generate next sibling at the
    // next parent cut node (next sibling of grandparent)

    // if expanding sibling of a non-leaf node, always expand the next sibling as CUT node... comparing it with
    // the current best in subtree

    printf("\nTotal nodes: %d\n", nCurr);
    printf("\nAfter scan, rejected: %d, to be expanded: %d\n", nRejected, nExpnded);

    free (fullCurrentFrontier);
    return 0.0;
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


int main()
{
    printf("\n\nSize of node is %d bytes\n\n", sizeof(Node));

    int depth = 6;

    //srand(time(NULL));
    // 3, and 4 are good. 6,7,8 are very good
    // 10 is excellent
    srand(10);
    // generate a random tree
    printf ("generating random tree of depth %d\n", depth);
    Node root = {0};
    START_TIMER
    genTree(&root, depth);
    STOP_TIMER
    printf("random tree generated, total nodes: %d, leaf nodes: %d, time: %g ms\n", gTotalNodes, gLeafNodes, gTime);



    float bestVal = 0;
    // search the best move using min-max search
    printf("searching the tree using min-max\n");
    START_TIMER
    bestVal = negaMax(&root, depth, depth);
    STOP_TIMER
    printf ("best move %d, score: %f, time: %g\n", root.bestChild, root.nodeVal, gTime);

    // search the best move using alpha-beta search
    printf("searching the tree using alpha-beta\n");
    START_TIMER
    bestVal = alphabeta(&root, depth, depth, -INF, INF);
    STOP_TIMER
    printf ("best move %d, score: %f\nnodes visited (leaves/interior/total): %d/%d/%d\n", 
            root.bestChild, root.nodeVal, gLeafNodesVisited, gInteriorNodesVisited, gLeafNodesVisited + gInteriorNodesVisited);
    printf("time taken: %g\n", gTime);

    
    START_TIMER
    exploreTree(&root, depth);
    STOP_TIMER
    printf("time taken: %g\n", gTime);    
    

    float val;
    START_TIMER
    val = SSS_star(&root, depth);
    STOP_TIMER
    printf("SSS* best node: %d, score: %f, nodes explored: %d, time taken: %g\n", root.bestChild, val, g_sssNodes, gTime);

    getchar();

    return 0;
}