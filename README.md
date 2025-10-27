# DFI-CPPI
A multi-omics framework integrating proteomic  correlations, partial correlations, and Gene Ontology semantic similarities to construct  context-specific PPI networks.
<img width="1685" height="1019" alt="image" src="https://github.com/user-attachments/assets/f9d10eef-d792-41ba-97df-ae1c26ad9475" />
calculate_cor_prob.ipynb
  Step 1: Data Loading and Preprocessing  
  Load proteomics data (CSV files).  
  Calculate correlation coefficients (Pearson and Partial Correlation).  
  Flatten the upper triangular part of the correlation coefficient matrices into a list of protein pairs.  

  Step 2: Logistic Regression Model Training and Prediction  
  Load known protein-protein interactions (from the CORUM database) as labels.  
  Train a Logistic Regression model.  
  Output prediction probabilities.  

  Step 3: XGBoost Model Training and Evaluation  
  Train XGBoost models using various feature combinations and perform hyperparameter tuning.  
  Evaluate the performance of different feature combinations (ROC curves, PR curves, etc.).  
  Save the best model and prediction results.  

  Step 4: Result Saving
calculate_cor_prob.ipynb
