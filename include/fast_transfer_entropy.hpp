/** MIT LICENSE
 * MIT License

 * Copyright (c) 2025 Corey Joshua Herr
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file fast_transfer_entropy.hpp
 * @brief A header only library for calculating transfer entropy between two symbol sequences.
 */

#ifndef PROBABILITY_TREE_H
#define PROBABILITY_TREE_H

#include <cmath>
#include <forward_list>
#include <iostream>
#include <vector>

namespace fast_entropy {

// FORWARD DECLARATIONS

/**
 * @brief Represents a node in the probability tree.
 */
struct Node
{
  int symbol; 
  int count = 0; 
  Node* parent = nullptr; 
  std::forward_list<Node*> children; 
};

using std::vector;
using Symbol = int;
using SymbolSequence = vector<Symbol>;

/**
 * @brief A class representing a probability tree for symbol sequences.
 */
class ProbabilityTree
{
public:
  /**
   * @brief Constructs a ProbabilityTree.
   * @param numSymbols Number of unique symbols.
   * @param numLevels Depth of the tree.
   * @param maxExpectedNodes Maximum number of expected nodes.
   */
  ProbabilityTree(int numSymbols, int numLevels, int maxExpectedNodes);

  /**
   * @brief Destructor to clean up allocated memory.
   */
  ~ProbabilityTree();

  /**
   * @brief Inserts a symbol sequence into the tree.
   * @param symbol The symbol sequence to insert.
   */
  void insertSymbol(SymbolSequence& symbol);

  /**
   * @brief Calculates the entropy of the tree.
   * @return The calculated entropy.
   */
  double calculateEntropy();

  /**
   * @brief Prints the entire tree structure.
   */
  void printTree();

  /**
   * @brief Recursively prints elements of the tree starting from a node.
   * @param Node The starting node.
   */
  void printElementRecursive(Node* Node);

  /**
   * @brief Finds the base count for a node at a specific level.
   * @param node The node to start from.
   * @param level The level to find the base count.
   * @return The base count.
   */
  int findBaseCount(Node* node, int level);

  /**
   * @brief Updates the base counts for all root nodes.
   */
  void updateBaseCount();

  // Total count of all symbols in the tree.
  int totalCount = 0;
  // List of leaf nodes in the tree.
  vector<Node*> leafNodes;
  // List of root nodes in the tree.
  vector<Node*> rootNodes;

  /**
   * @brief Calculates the KL Transfer Entropy between two trees.
   * @param treeXY The first tree.
   * @param treeYX The second tree.
   * @param delta The time delay.
   * @return The calculated KL Transfer Entropy.
   */
  friend double calculateSymbolicTransferEntropy(ProbabilityTree& treeXY,
                                           ProbabilityTree& treeYX,
                                           const int delta);

private:
  /**
   * @brief Recursively adds a symbol sequence to the tree with checks.
   * @param symbol The symbol sequence to add.
   * @param level The current level in the tree.
   * @param currentNode The current node being processed.
   */
  void addRecursiveCheck(SymbolSequence& symbol, int level, Node* currentNode);

  /**
   * @brief Adds a symbol sequence to the tree without checks.
   * @param symbol The symbol sequence to add.
   * @param level The current level in the tree.
   * @param currentNode The current node being processed.
   */
  void addNoCheck(SymbolSequence& symbol, int level, Node* currentNode);

  // 2D vector representing the tree structure.
  vector<vector<Node*>> tree;
  // Number of unique symbols.
  int numSymbols_;
  // Depth of the tree.
  int numLevels_;
  // Maximum number of expected nodes.
  int maxExpectedNodes_;
};

/**
 * @brief Calculates the KL Transfer Entropy between two symbol sequences.
 * @param X The first symbol sequence.
 * @param Y The second symbol sequence.
 * @param delta The time delay.
 * @param numSymbols The number of unique symbols.
 * @return The calculated KL Transfer Entropy.
 */
double calculateSymbolicTransferEntropy(const SymbolSequence& X,
                                  const SymbolSequence& Y,
                                  const int delta,
                                  const int numSymbols);

//
// IMPLEMENTATIONS BELOW
//

/**
 * @brief Finds the ancestor of a node at a specific depth.
 * @param node The starting node.
 * @param k The depth to find the ancestor.
 * @return The ancestor node.
 */
inline Node* findAncestor(Node* node, int k)
{
  Node* parentNode = node;
  for (int i = 0; i < k; i++) {
    parentNode = parentNode->parent;
  }
  return parentNode;
}

/**
 * @brief Calculates the KL Transfer Entropy between two probability trees.
 * @param treeXY The first tree.
 * @param treeYX The second tree.
 * @param delta The time delay.
 * @return The calculated KL Transfer Entropy.
 */
inline double calculateSymbolicTransferEntropy(ProbabilityTree& treeXY,
                                         ProbabilityTree& treeYX,
                                         const int delta)
{
  // Calculate the Transfer entropy given two probability trees
  // sum(p(yn,yk,xk)*log2(p(yn|yk,xk)/p(yn|yk)))
  // HXY calculates sum(p(yn,yk,xk)*log2(p(yn|yk,xk))
  // HYX calculates sum(p(yn,yk,xk)*log2(p(yn|yk))

  double H = 0;

  // Loop to calcualte HXY
  for (Node* node : treeXY.leafNodes) {
    double currentCount = node->count;
    double parentCount = node->parent->count;
    H += (currentCount / treeXY.totalCount) * log2(currentCount / parentCount);
  }

  // Loop to calculate HYX
  // This loop is pretty expensive because of findAncestor
  // TODO:See if there is any way to optimize data structure for these calls
  for (Node* node : treeYX.leafNodes) {
    double currentCount = node->count;
    Node* ancestorN = findAncestor(node, delta - 1);
    double countYn = ancestorN->count;
    double countYnk = ancestorN->parent->count;
    H -= (currentCount / treeYX.totalCount) * log2(countYn / countYnk);
  }
  return H;
}

/**
 * @brief Calculates the Transfer Entropy between two probability trees.
 * @param treeY The first tree.
 * @param treeXY The second tree.
 * @param delta The time delay.
 * @return The calculated Transfer Entropy.
 */

/**
 * @brief Calculates the KL Transfer Entropy between two symbol sequences.
 * @param X The first symbol sequence.
 * @param Y The second symbol sequence.
 * @param delta The time delay.
 * @param numSymbols The number of unique symbols.
 * @return The calculated KL Transfer Entropy.
 */
inline double calculateSymbolicTransferEntropy(const SymbolSequence& X,
                                         const SymbolSequence& Y,
                                         const int delta,
                                         const int numSymbols)
{
  int T = X.size();
  double TXY = 0;
  double TYX = 0;

  {
    ProbabilityTree treeXY = ProbabilityTree(numSymbols, 2 * delta - 1, T);
    ProbabilityTree treeYX = ProbabilityTree(numSymbols, 2 * delta - 1, T);
    for (int i = delta; i < X.size(); i++) {
      vector<int> symbolXY(X.begin() + i - delta, X.begin() + i - 1);
      symbolXY.insert(symbolXY.end(), Y.begin() + i - delta, Y.begin() + i);
      [[assume(symbolXY.size() == 2 * delta - 1)]];

      vector<int> symbolYX(Y.begin() + i - delta, Y.begin() + i);
      symbolYX.insert(symbolYX.end(), X.begin() + i - delta, X.begin() + i - 1);
      [[assume(symbolYX.size() == 2 * delta - 1)]];

      treeXY.insertSymbol(symbolXY);
      treeYX.insertSymbol(symbolYX);
    }
    if (delta == 2) {
      treeXY.updateBaseCount();
      treeYX.updateBaseCount();
    }

    TXY = calculateSymbolicTransferEntropy(treeXY, treeYX, delta);
  }
  {
    ProbabilityTree treeYX(numSymbols, 2 * delta - 1, T);
    ProbabilityTree treeXY(numSymbols, 2 * delta - 1, T);
    for (int i = delta; i < T; i++) {
      vector<int> symbolYX(Y.begin() + i - delta, Y.begin() + i - 1);
      symbolYX.insert(symbolYX.end(), X.begin() + i - delta, X.begin() + i);

      vector<int> symbolXY(X.begin() + i - delta, X.begin() + i);
      symbolXY.insert(symbolXY.end(), Y.begin() + i - delta, Y.begin() + i - 1);

      treeYX.insertSymbol(symbolYX);
      treeXY.insertSymbol(symbolXY);
    }
    if (delta == 2) {
      treeXY.updateBaseCount();
      treeYX.updateBaseCount();
    }

    TYX = calculateSymbolicTransferEntropy(treeYX, treeXY, delta);
  }

  return TXY - TYX;
}


/**
 * @brief Deletes all nodes in the tree using DFS.
 * @param node The starting node.
 */
inline void DFSDelete(Node* node)
{
  for (auto& child : node->children) {
    if (child)
      DFSDelete(child);
  }
  delete node;
  node = nullptr;
}

/**
 * @brief Destructor for the ProbabilityTree class.
 */
inline ProbabilityTree::~ProbabilityTree()
{
  for (int i = 0; i < numSymbols_; i++) {
    for (int j = 0; j < numSymbols_; j++) {
      if (tree[i][j]) {
        DFSDelete(tree[i][j]);
      }
    }
  }
  for (int i = 0; i < numSymbols_; i++) {
    if (rootNodes[i]) {
      delete rootNodes[i];
      rootNodes[i] = nullptr;
    }
  }
}

/**
 * @brief Prints the entire tree structure.
 */
inline void ProbabilityTree::printTree()
{
  for (int i = 0; i < numSymbols_; i++) {
    for (int j = 0; j < numSymbols_; j++) {
      if (tree[i][j]) {
        printElementRecursive(tree[i][j]);
      }
    }
  }
}

/**
 * @brief Recursively prints elements of the tree starting from a node.
 * @param Node The starting node.
 */
inline void ProbabilityTree::printElementRecursive(Node* Node)
{
  for (auto const& child : Node->children) {
    printElementRecursive(child);
  }
  std::cout << "Node: " << Node->symbol << " count: " << Node->count << "\n";
}

/**
 * @brief Constructs a ProbabilityTree.
 * @param numSymbols Number of unique symbols.
 * @param numLevels Depth of the tree.
 * @param maxExpectedNodes Maximum number of expected nodes.
 */
inline ProbabilityTree::ProbabilityTree(int numSymbols,
                                        int numLevels,
                                        int maxExpectedNodes)
{
  numSymbols_ = numSymbols;
  numLevels_ = numLevels;
  maxExpectedNodes_ = maxExpectedNodes;
  tree =
    vector<vector<Node*>>(numSymbols_, vector<Node*>(numSymbols_, nullptr));
  rootNodes = vector<Node*>(numSymbols_, nullptr);
  for (int i = 0; i < numSymbols_; i++) {
    rootNodes[i] = new Node;
    rootNodes[i]->symbol = i;
    rootNodes[i]->count = 0;
  }
  leafNodes.reserve(maxExpectedNodes_);
}

/**
 * @brief Adds a symbol sequence to the tree without checks.
 * @param symbol The symbol sequence to add.
 * @param level The current level in the tree.
 * @param currentNode The current node being processed.
 */
inline void ProbabilityTree::addNoCheck(SymbolSequence& symbol,
                                        int level,
                                        Node* currentNode)
{
  Node* node = currentNode;
  for (int i = level; i < numLevels_; i++) {
    Node* newChildNode = new Node{ symbol[i], 1 };
    node->children.push_front(newChildNode);

    if (i == numLevels_ - 1)
      leafNodes.push_back(newChildNode);

    newChildNode->parent = node;

    node = newChildNode;
  }
}

/**
 * @brief Recursively adds a symbol sequence to the tree with checks.
 * @param symbol The symbol sequence to add.
 * @param level The current level in the tree.
 * @param currentNode The current node being processed.
 */
inline void ProbabilityTree::addRecursiveCheck(SymbolSequence& symbol,
                                               int level,
                                               Node* currentNode)
{
  if (level == numLevels_) {
    return;
  }
  for (auto& child : currentNode->children) {
    if (child->symbol == symbol[level]) {
      child->count++;
      addRecursiveCheck(symbol, level + 1, child);
      return;
    }
  }
  addNoCheck(symbol, level, currentNode);
}

/**
 * @brief Calculates the entropy of the tree.
 * @return The calculated entropy.
 */
inline double ProbabilityTree::calculateEntropy()
{
  updateBaseCount();
  double entropy = 0;
  for (int i = 0; i < numSymbols_; i++) {
    if (rootNodes[i]->count != 0) {
      double p = rootNodes[i]->count / (double)totalCount;
      entropy -= p * log2(p);
    }
  }
  return entropy;
}

/**
 * @brief Inserts a symbol sequence into the tree.
 * @param symbol The symbol sequence to insert.
 */
inline void ProbabilityTree::insertSymbol(SymbolSequence& symbol)
{
  // symbol must be length >2
  totalCount++;
  if (!tree[symbol[0]][symbol[1]]) {
    tree[symbol[0]][symbol[1]] = new Node;
    tree[symbol[0]][symbol[1]]->symbol = symbol[0];
    tree[symbol[0]][symbol[1]]->count = 0;
    tree[symbol[0]][symbol[1]]->parent = rootNodes[symbol[0]];
  }
  Node* currentNode = tree[symbol[0]][symbol[1]];
  currentNode->count++;
  addRecursiveCheck(symbol, 2, currentNode);
}

/**
 * @brief Updates the base counts for all root nodes.
 */
inline void ProbabilityTree::updateBaseCount()
{
  for (int i = 0; i < numSymbols_; i++) {
    rootNodes[i]->count = 0;
  }

  for (int i = 0; i < numSymbols_; i++) {
    for (int j = 0; j < numSymbols_; j++) {
      Node* child = tree[i][j];
      if (child) {
        child->parent->count += child->count;
      }
    }
  }
}

} // namespace fast_entropy

#endif // !PROBABILITY_TREE_H
