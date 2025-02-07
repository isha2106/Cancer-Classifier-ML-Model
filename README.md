Detection of Cancer Type based on Gene Expression Profiles

## Introduction
In bioinformatics, there is significant use of gene expression profiles for numerous number of studies. From
identification of cancer to finding therapeutic targets gene expression data has proven to be of utmost
importance. The primary goal of this project is to use the gene expression data for categorizing different
types of cancer. Here, machine learning techniques are used for the same which can be eﬀicient and a faster
approach for cancer type identification. This project could pave the way for more personalized medicine
approaches and also contribute to ongoing efforts in oncology research and treatment strategies.

## Objective
To determine the type of cancer based on the gene expression profile using machine learning
models.

## Formulated data mining approach
This is a supervised learning problem where the input variables are the gene expression data and the target
variable is the type of cancer. The challenge involves processing the high dimensional dataset, selecting
the most relevant features and choosing machine learning algorithms that can effectively classify the cancer
types. Additionally, the project will explore the potential of combining models through ensemble methods
to enhance predictive accuracy and robustness.

## Requirements
RNA-Seq (HiSeq) PANCAN dataset with instances (patients) and features (gene expressions).

## Outline
This project uses the RNA-Seq (HiSeq) PANCAN dataset from the UC Irvine Machine Learning Repository
to determine the type of cancer (BRCA, KIRC, COAD, LUAD, PRAD) using gene expressions from patients.
The classification model uses supervised machine learning approach. The project will explore Support Vector
Machine (SVM), Random Forest, and Gradient Boosting algorithms known for their high-dimensional data
handling, feature interaction capture, and skewness mitigation capabilities. Utilizing data normalization
techniques and feature selection, the project aims for precise cancer classification, evaluated using cross-
validation and accuracy/F1 scores.
