#include <iostream>
#include <vector>
#include <cmath>
#include<algorithm>

using namespace std;

struct Node {
    int featureIndex;
    double threshold;
    bool isLeaf;
    bool decision;
    Node* left;
    Node* right;
};

Node* head;
/*struct Node *create_newnode(int featureindex,double Threshold,bool isleaf,bool Decision)
{
    struct Node *newnode = (struct Node *)malloc(sizeof(struct Node));
    newnode->featureIndex = featureindex;
    newnode->left = newnode->right = NULL;
    return newnode;
}*/

double calculateEntropy(vector<bool>& labels) {
    int totalInstances = labels.size();
    int trueCount = 0;

    for (bool label : labels) {
        if (label) {
            trueCount++;
        }
    }

    double trueProbability = (double)(trueCount) / totalInstances;
    double falseProbability = 1.0 - trueProbability;

    // Avoid log(0)
    if (trueProbability == 0.0 || falseProbability == 0.0) {
        return 0.0;
    }

    return -trueProbability * log2(trueProbability) - falseProbability * log2(falseProbability);
}

bool areAllLabelsEqual(vector<bool>& labels) {
    if (labels.empty()) {
        return true;  // Empty vector is considered to have equal labels
    }

    bool firstLabel = labels[0];
    for (bool label : labels) {
        if (label != firstLabel) {
            return false;  // Labels are not equal
        }
    }

    return true;  // All labels are equal
}

pair<int, double> findBestSplit(vector<vector<double>>& data, vector<bool>& labels) {
    int numFeatures = data[0].size() - 1;  // Last column is the label
    int numInstances = data.size();

    double currentEntropy = calculateEntropy(labels);
    double maxInfoGain = 0.0;
    int bestFeatureIndex = -1;
    double bestThreshold = 0.0;

    for (int featureIndex = 0; featureIndex < numFeatures; ++featureIndex) {
        // Sort data based on the current feature
        vector<pair<double, bool>> featureData;
        for (int instance = 0; instance < numInstances; ++instance) {
            featureData.push_back({data[instance][featureIndex], labels[instance]});
        }

        sort(featureData.begin(), featureData.end());

        // Iterate through possible thresholds
        for (int i = 0; i < numInstances - 1; ++i) {

            double threshold = (featureData[i].first + featureData[i + 1].first) / 2.0;

            vector<bool> leftLabels(featureData.size(), false);
            vector<bool> rightLabels(featureData.size(), false);

            // Update left and right labels based on threshold
            for (int j = 0; j <= i; ++j) {
                leftLabels[j] = featureData[j].second;
            }

            for (int j = i + 1; j < numInstances; ++j) {
                rightLabels[j] = featureData[j].second;
            }

            double infoGain = currentEntropy - (calculateEntropy(leftLabels) * (i + 1) +
                                                 calculateEntropy(rightLabels) * (numInstances - i - 1)) / numInstances;

            if (infoGain > maxInfoGain) {
                maxInfoGain = infoGain;
                bestFeatureIndex = featureIndex;
                bestThreshold = threshold;
            }
        }
    }

    return {bestFeatureIndex, bestThreshold};
}

void deleteTree(Node* root) {
    if (root == nullptr) {
        return;
    }

    deleteTree(root->left);
    deleteTree(root->right);
    delete root;
}

Node* buildDecisionTree(vector<vector<double>>& data, vector<bool>& labels) {
    Node* root = new Node;

    for(bool it : labels)
    {
      cout << it << " ";
    }cout << "\n";

    // Base case: If all labels are the same, create a leaf node
    if (areAllLabelsEqual(labels)) {
        root->isLeaf = true;
        root->decision = labels[0];
        cout << "areAllLabelsfunction called\n";
        return root;
    }

    // Find the best split
    pair<int, double> bestSplit = findBestSplit(data, labels);
    root->featureIndex = bestSplit.first;
    root->threshold = bestSplit.second;

    // Split the data into left and right subsets
    vector<std::vector<double>> leftData, rightData;
    vector<bool> leftLabels, rightLabels;

    for (int i = 0; i < data.size(); ++i) {
        if (data[i][bestSplit.first] <= bestSplit.second) {
            leftData.push_back(data[i]);
            leftLabels.push_back(labels[i]);
        } else {
            rightData.push_back(data[i]);
            rightLabels.push_back(labels[i]);
        }
    }

    // Recursively build left and right subtrees
    root->left = buildDecisionTree(leftData, leftLabels);
    root->right = buildDecisionTree(rightData, rightLabels);

    return root;
}

bool predict(Node* root, vector<double>& sample) {
    if (root->isLeaf) {
        return root->decision;
    }

    if (sample[root->featureIndex] <= root->threshold) {
        return predict(root->left, sample);
    } else {
        return predict(root->right, sample);
    }
}

void printTree(Node *head)
{
  if(head->isLeaf)
  {
     cout << "Leaf reached\n";
     cout <<  head->decision << "\n" << head->featureIndex << "\n" << head->isLeaf <<"\n";
     return;
  }
  cout << "Tree node info:\n";
  cout <<  head->featureIndex << "\n" << head->threshold << "\n" << head->decision<<"\n";
  printTree(head->left);
  printTree(head->right);

}
int main() {
    // Demo dataset
    vector<vector<double>> data = {{1.1, 1, 1, 1.1},
                                   {1, 0.1, 1, 2.0},
                                   {1.5, 1.1, 1, 1.0}};
    vector<bool> labels = {true, true, false};

    // Build decision tree
    Node* decisionTree = buildDecisionTree(data, labels);

    head = decisionTree;

    printTree(head);

    // Make predictions
   // vector<double> sample = {1.2, 0.9, 1, 1.2};  // Replace with your own sample
   // bool prediction = predict(decisionTree, sample);

    // Clean up memory
    //deleteTree(decisionTree);

    //cout << decisionTree->threshold << "\n"<< decisionTree->left->threshold;

    //cout << "Prediction: " << prediction << std::endl;

    return 0;
}
