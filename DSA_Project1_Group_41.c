#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define maxchild 4
#define minchild 2

typedef struct
{
    double x;
    double y;
} Point;

typedef struct
{
    double xmax;
    double xmin;
    double ymax;
    double ymin;
} Mbr;

typedef struct node
{
    int isleaf;
    int count;           // number of children or points
    struct node *parent; // address of parent node
    struct node **child; // array of children addresses
    Mbr rectangle;       // mbr of node
    Point **data;        // array of points
} Node;

typedef struct
{
    Node *root;
} RTree;

Node *createNode(int isleaf)
{
    Node *myNode = (Node *)malloc(sizeof(Node));
    myNode->count = 0;
    myNode->parent = NULL;
    myNode->child = (Node **)malloc((maxchild + 1) * sizeof(Node *));
    myNode->isleaf = isleaf;
    for (int i = 0; i < maxchild + 1; i++)
    {
        myNode->child[i] = NULL;
    }
    myNode->rectangle.xmax = 0.0;
    myNode->rectangle.xmin = 0.0;
    myNode->rectangle.ymax = 0.0;
    myNode->rectangle.ymin = 0.0;
    myNode->data = (Point **)malloc((maxchild + 1) * sizeof(Point *));
    for (int i = 0; i < maxchild + 1; i++)
    {
        myNode->data[i] = NULL;
    }
    return myNode;
}

int contains_point(Mbr mbr, Point point)
{
    return (mbr.xmin <= point.x && mbr.xmax >= point.x && mbr.ymin <= point.y && mbr.ymax >= point.y);
}

//  calulate mbr of a node if point is added
Mbr getNewMbr(Mbr mbr, Point *point)
{
    Mbr newMbr;
    newMbr.xmin = fmin(mbr.xmin, point->x);
    newMbr.xmax = fmax(mbr.xmax, point->x);
    newMbr.ymin = fmin(mbr.ymin, point->y);
    newMbr.ymax = fmax(mbr.ymax, point->y);
    return newMbr;
}

// calculate mbr of a node if mbr is added
Mbr getNewMbr2(Mbr mbr, Mbr mbr2)
{
    Mbr newMbr;
    newMbr.xmin = fmin(mbr.xmin, mbr2.xmin);
    newMbr.xmax = fmax(mbr.xmax, mbr2.xmax);
    newMbr.ymin = fmin(mbr.ymin, mbr2.ymin);
    newMbr.ymax = fmax(mbr.ymax, mbr2.ymax);
    return newMbr;
}
// update mbr using the child nodes
Mbr updateMbr(Node *node)
{
    Mbr mbr;
    mbr.xmin = node->child[0]->rectangle.xmin;
    mbr.xmax = node->child[0]->rectangle.xmax;
    mbr.ymin = node->child[0]->rectangle.ymin;
    mbr.ymax = node->child[0]->rectangle.ymax;
    for (int i = 1; i < node->count; i++)
    {
        if (node->child[i] == NULL)
            continue;
        mbr.xmin = fmin(mbr.xmin, node->child[i]->rectangle.xmin);
        mbr.xmax = fmax(mbr.xmax, node->child[i]->rectangle.xmax);
        mbr.ymin = fmin(mbr.ymin, node->child[i]->rectangle.ymin);
        mbr.ymax = fmax(mbr.ymax, node->child[i]->rectangle.ymax);
    }
    return mbr;
}
//  calculate area of mbr
double calculateArea(Mbr mbr)
{
    return (mbr.xmax - mbr.xmin) * (mbr.ymax - mbr.ymin);
}

// search for a point in the tree
void searchRtree(Node *node, Point *searchPoint)
{
    if (node->isleaf)
    {
        // Search leaf node
        for (int i = 0; i < node->count; i++)
        {
            if (node->data[i] == NULL)
                continue;
            if (node->data[i]->x == searchPoint->x && node->data[i]->y == searchPoint->y)
            {
                printf("Found point (%f, %f)", node->data[i]->x, node->data[i]->y);
            }
        }
    }
    else
    {
        // Search non-leaf node
        for (int i = 0; i < node->count; i++)
        {
            if (contains_point(node->child[i]->rectangle, *searchPoint))
            {
                searchRtree(node->child[i], searchPoint);
            }
        }
    }
}

// pick the 2 children with most inefficient area
void pickSeeds(Node *node, int *indexOfSeed1, int *indexOfSeed2)
{
    double maxInefficiency = 0.0;
    double d;
    for (int i = 0; i < node->count; i++)
    {
        for (int j = i + 1; j < node->count; j++)
        {
            Mbr mbr;
            mbr.xmin = fmin(node->child[i]->rectangle.xmin, node->child[j]->rectangle.xmin);
            mbr.xmax = fmax(node->child[i]->rectangle.xmax, node->child[j]->rectangle.xmax);
            mbr.ymin = fmin(node->child[i]->rectangle.ymin, node->child[j]->rectangle.ymin);
            mbr.ymax = fmax(node->child[i]->rectangle.ymax, node->child[j]->rectangle.ymax);
            d = calculateArea(mbr) - calculateArea(node->child[i]->rectangle) - calculateArea(node->child[j]->rectangle);
            if (d > maxInefficiency)
            {
                *indexOfSeed1 = i;
                *indexOfSeed2 = j;
                maxInefficiency = d;
            }
        }
    }
}

// pick the children
void pickNext(Node *node, Node *newNode1, Node *newNode2)
{
    if (node->isleaf)
    {
        for (int i = 0; i < node->count; i++)
        {
            if (node->data[i] != NULL)
            {
                Mbr new_mbr1 = getNewMbr(newNode1->rectangle, node->data[i]);
                Mbr new_mbr2 = getNewMbr(newNode2->rectangle, node->data[i]);
                double area1 = calculateArea(new_mbr1) - calculateArea(newNode1->rectangle);
                double area2 = calculateArea(new_mbr2) - calculateArea(newNode2->rectangle);
                if (area1 < area2 || (area1 == area2 && calculateArea(newNode1->rectangle) < calculateArea(newNode2->rectangle)))
                {
                    newNode1->data[newNode1->count] = node->data[i];
                    newNode1->rectangle = new_mbr1;
                    newNode1->count++;
                }
                else
                {
                    newNode2->data[newNode2->count] = node->data[i];
                    newNode2->rectangle = new_mbr2;
                    newNode2->count++;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < node->count; i++)
        {
            if (node->child[i] == NULL)
                continue;
            Mbr mbr1 = getNewMbr2(newNode1->rectangle, node->child[i]->rectangle);
            Mbr mbr2 = getNewMbr2(newNode2->rectangle, node->child[i]->rectangle);
            double area_inc1 = calculateArea(mbr1) - calculateArea(newNode1->rectangle);
            double area_inc2 = calculateArea(mbr2) - calculateArea(newNode2->rectangle);
            if (area_inc1 < area_inc2 || (area_inc1 == area_inc2 && newNode1->count < newNode2->count))
            {
                newNode1->child[newNode1->count] = node->child[i];
                newNode1->rectangle = updateMbr(newNode1);
                newNode1->count++;
            }
            else
            {
                newNode2->child[newNode2->count] = node->child[i];
                newNode2->rectangle = updateMbr(newNode2);
                newNode2->count++;
            }
        }
    }
}

void adjustTree(Node *node)
{
    // checks to see that the node is not a root
    if (node->parent != NULL)
    {
        node->parent->rectangle = updateMbr(node->parent); // adjusts covering entry in parent MBR
        adjustTree(node->parent);                          // moves up to the next level
    }
}

// split the node into two using quadratic split
void splitNode(Node *node, RTree *rtree)
{
    Node *newNode1 = createNode(node->isleaf);
    Node *newNode2 = createNode(node->isleaf);
    if (node->isleaf)
    {
        // find two points that are farthest apart
        double eulerDistance, maxDistance = -1.0;
        int indexOfSeed1, indexOfSeed2;
        for (int i = 0; i < node->count; i++)
        {
            for (int j = i + 1; j < node->count; j++)
            {
                if (node->data[i] == NULL || node->data[j] == NULL)
                    continue;
                eulerDistance = sqrt(pow(node->data[i]->x - node->data[j]->x, 2) + pow(node->data[i]->y - node->data[j]->y, 2));
                if (eulerDistance > maxDistance)
                {
                    maxDistance = eulerDistance;
                    indexOfSeed1 = i;
                    indexOfSeed2 = j;
                }
            }
        }
        // add farthest points to new nodes
        newNode1->data[newNode1->count] = node->data[indexOfSeed1];
        newNode1->count++;
        newNode2->data[newNode2->count] = node->data[indexOfSeed2];
        newNode2->count++;

        // remove farthest points from node
        node->data[indexOfSeed1] = NULL;
        node->data[indexOfSeed2] = NULL;

        // calculate mbr of new nodes
        newNode1->rectangle.xmax = newNode1->data[0]->x;
        newNode1->rectangle.xmin = newNode1->data[0]->x;
        newNode1->rectangle.ymax = newNode1->data[0]->y;
        newNode1->rectangle.ymin = newNode1->data[0]->y;
        newNode2->rectangle.xmax = newNode2->data[0]->x;
        newNode2->rectangle.xmin = newNode2->data[0]->x;
        newNode2->rectangle.ymax = newNode2->data[0]->y;
        newNode2->rectangle.ymin = newNode2->data[0]->y;

        // add other points to newNode1 or newNode2
        pickNext(node, newNode1, newNode2);
        // set the child pointer of parent node to new child nodes
        if (node->parent == NULL)
        {
            Node *newRoot = createNode(0);
            newRoot->child[0] = newNode1;
            newRoot->child[1] = newNode2;
            newRoot->rectangle = updateMbr(newRoot);
            newRoot->count = 2;
            newNode1->parent = newRoot;
            newNode2->parent = newRoot;
            rtree->root = newRoot;
        }
        else
        {
            // remove the node from parent and shift all nodes to remove gap in array
            for (int i = 0; i < node->parent->count; i++)
            {
                if (node->parent->child[i] == node)
                {
                    for (int j = i; j < node->parent->count - 1; j++)
                    {
                        node->parent->child[j] = node->parent->child[j + 1];
                    }
                    node->parent->child[node->parent->count - 1] = NULL;
                    node->parent->count--;
                    break;
                }
            }
            // if the count of parent is less that 4 then add the new nodes to parent where the child array has NULL
            if (node->parent->count < maxchild)
            {
                node->parent->child[node->parent->count] = newNode1;
                node->parent->child[node->parent->count + 1] = newNode2;
                node->parent->count = node->parent->count + 2;
                node->parent->rectangle = updateMbr(node->parent);
                newNode1->parent = node->parent;
                newNode2->parent = node->parent;
            }
            else
            {
                // if the count of parent is 4 then split the parent
                node->parent->child[node->parent->count] = newNode1;
                node->parent->child[node->parent->count + 1] = newNode2;
                node->parent->count = node->parent->count + 2;
                splitNode(node->parent, rtree);
            }
        }
    }
    else
    {
        int indexOfSeed1, indexOfSeed2;
        pickSeeds(node, &indexOfSeed1, &indexOfSeed2);
        newNode1->child[newNode1->count] = node->child[indexOfSeed1];
        newNode1->count++;
        newNode2->child[newNode2->count] = node->child[indexOfSeed2];
        newNode2->count++;
        newNode1->rectangle = updateMbr(newNode1);
        newNode2->rectangle = updateMbr(newNode2);

        // remove both seeds from node and shift the child array to remove gap
        for (int i = 0; i < node->count; i++)
        {
            if (node->child[i] == node->child[indexOfSeed1] || node->child[i] == node->child[indexOfSeed2])
            {
                for (int j = i; j < node->count - 1; j++)
                {
                    node->child[j] = node->child[j + 1];
                }
                node->child[node->count - 1] = NULL;
                node->count--;
            }
        }

        // distribute remaining childre to newNode1 and newNode2 based on minimum area increase on adding new child
        pickNext(node, newNode1, newNode2);
        // set the child pointer of parent node to new child nodes
        if (node->parent == NULL)
        {
            Node *newRoot = createNode(0);
            newRoot->child[0] = newNode1;
            newRoot->child[1] = newNode2;
            newRoot->rectangle = updateMbr(newRoot);
            newRoot->count = 2;
            newNode1->parent = newRoot;
            newNode2->parent = newRoot;
            rtree->root = newRoot;
        }
        else
        {
            // remove the node from parent and shift all nodes to remove gap in array
            for (int i = 0; i < node->parent->count; i++)
            {
                if (node->parent->child[i] == node)
                {
                    for (int j = i; j < node->parent->count - 1; j++)
                    {
                        node->parent->child[j] = node->parent->child[j + 1];
                    }
                    node->parent->child[node->parent->count - 1] = NULL;
                    node->parent->count--;
                    break;
                }
            }
            // if the count of parent is less that 4 then add the new nodes to parent where the child array has NULL
            if (node->parent->count < maxchild)
            {
                node->parent->child[node->parent->count] = newNode1;
                node->parent->child[node->parent->count + 1] = newNode2;
                node->parent->count = node->parent->count + 2;
                node->parent->rectangle = updateMbr(node->parent);
                newNode1->parent = node->parent;
                newNode2->parent = node->parent;
            }
            else
            {
                // if the count of parent is 4 then split the parent
                node->parent->child[node->parent->count] = newNode1;
                node->parent->child[node->parent->count + 1] = newNode2;
                node->parent->count = node->parent->count + 2;
                splitNode(node->parent, rtree);
            }
        }
    }
}

int areaIncrease(Node *node, Point *point)
{
    double increaseInArea[5];
    for (int i = 0; i < node->count; i++)
    {
        if (node->child[i] != NULL)
        {
            increaseInArea[i] = calculateArea(getNewMbr(node->child[i]->rectangle, point)) - calculateArea(node->child[i]->rectangle);
        }
    }
    // Choose the child node with the minimum area increase
    int j = 0;
    for (int i = 1; i < node->count; i++)
    {
        if (increaseInArea[i] < increaseInArea[j])
        {
            j = i;
        }
    }
    return j;
}

// choose the leaf node where the point should be inserted
Node *chooseLeaf(Node *node, Point *point)
{
    // If the node is a leaf node, return the node
    if (node->isleaf)
    {
        return node;
    }
    // If the node is not a leaf node, choose the appropriate child node
    else
    {
        // Calculate the area increase of each child node
        int minimumIndex = areaIncrease(node, point);
        return chooseLeaf(node->child[minimumIndex], point);
    }
}

// insert a point into the node
void insertData(Node *node, Point *point, RTree *rtree)
{
    // Check if the node is a leaf node
    if (node->isleaf)
    {
        // If the node is not full, insert the point into the node
        if (node->count < 4)
        {
            node->data[node->count] = point;
            node->rectangle = getNewMbr(node->rectangle, point);
            node->count++;
        }
        // If the node is full, split the node
        else
        {
            node->data[node->count] = point;
            node->count++;
            splitNode(node, rtree);
            adjustTree(node);
        }
    }
    // If the node is not a leaf node, insert the point into the appropriate child node
    else
    {
        Node *childNode = chooseLeaf(node, point); // Find position for new record
        insertData(childNode, point, rtree);       // adding point to node
        adjustTree(childNode);                     // propogating changes upward
    }
}

// create a rtree
RTree *createRTree(Point *points, int n)
{
    // Create a root node
    Node *root = createNode(1);
    // Create a RTree
    RTree *tree = (RTree *)malloc(sizeof(RTree));
    tree->root = root;
    // Insert all the points into the RTree
    for (int i = 0; i < n; i++)
    {
        insertData(tree->root, &points[i], tree);
    }
    return tree;
}

void print_leaf(Node *node, int i)
{
    if (node->data[i] != NULL)
        printf("Leaf node with point (%f, %f)\n", node->data[i]->x, node->data[i]->y); // print points
}

// preorder traversal of the tree
void preorderTraversal(Node *node)
{
    static int flag = 0;
    if (node != NULL)
    {
        if (flag == 0)
        {
            printf("Root Node with %d children\n", node->count);
            flag++;
        }
        else
            printf("Internal node with %d children\n", node->count);
        for (int i = 0; i < node->count; i++)
        {
            if (node->isleaf)
            {
                print_leaf(node, i);
            }
            else
            {
                if (node->child[i] != NULL)
                {
                    printf("Internal node with mbr (xmin=%f, xmax=%f, ymin=%f, ymax=%f)\n", node->child[i]->rectangle.xmin, node->child[i]->rectangle.xmax, node->child[i]->rectangle.ymin, node->child[i]->rectangle.ymax);
                    preorderTraversal(node->child[i]); // recursive call
                }
            }
        }
    }
}

// main function which would read from a file and create a rtree
int main()
{
    // read from data.txt file
    FILE *fp;
    char *filename = (char *)malloc(sizeof(char) * 100);
    printf("Enter the filename: ");
    scanf("%s", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("File not found\n");
        return 0;
    }
    int n;
    printf("Enter the number of points: ");
    scanf("%d", &n);
    Point *points = (Point *)malloc(sizeof(Point) * n);
    for (int i = 0; i < n; i++)
    {
        fscanf(fp, "%lf %lf", &points[i].x, &points[i].y);
    }
    fclose(fp);
    // create a rtree
    RTree *tree = createRTree(points, n);
    // print the rtree
    Point p;
    preorderTraversal(tree->root);
    searchRtree(tree->root, &p);
    free(tree);
    free(points);
    return 0;
}
