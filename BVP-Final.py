#!/usr/bin/env python
# coding: utf-8

# In[28]:


from Bio import SeqIO
import pandas as pd
import subprocess
import os
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report


# In[2]:


sequences = []
labels = []

for record in SeqIO.parse("/Users/anandnambigovindarajan/Desktop/BVP/VFDB_setA_pro(core dataset).fas", "fasta"):  
    sequences.append(str(record.seq))
    labels.append(1) # 1 stands for virulent

# save to csv file
df.to_csv("virulence_dataset.csv", index=False)


# In[3]:


# Read in your negative sequence data
sequences_n = []
labels_n = []

for record in SeqIO.parse("/Users/anandnambigovindarajan/Desktop/BVP/negative dataset.fasta", "fasta"):  
    sequences_n.append(str(record.seq))
    labels_n.append(0) # 0 stands for non-virulent


# In[4]:


# Check the uniqueness of your sequences
unique_sequences = len(set(sequences))
print(f"Number of unique sequences: {unique_sequences} out of {len(sequences)} total sequences.")


# In[5]:


# Check the uniqueness of your sequences
unique_sequences_n = len(set(sequences_n))
print(f"Number of unique sequences: {unique_sequences_n} out of {len(sequences_n)} total sequences.")


# In[8]:


# Create a dataframe for the positive data
df = pd.DataFrame(list(zip(sequences, labels)), columns=["Sequence", "Label"])


# In[9]:


# Create a dataframe for the negative data
df_n = pd.DataFrame(list(zip(sequences_n, labels_n)), columns=["Sequence", "Label"])


# In[10]:


# Check the balance of your classes
class_counts = df['Label'].value_counts()
print(class_counts)


# In[11]:


# Check the balance of your classes
class_counts = df_n['Label'].value_counts()
print(class_counts)


# In[12]:


# Change permission for the iFeature.py
subprocess.call(["chmod", "+x", "/Users/anandnambigovindarajan/iFeature/iFeature.py"])


# In[17]:


#Add feature
df['Features'] = pd.Series([[] for _ in range(len(df))], index=df.index)
df_n['Features'] = pd.Series([[] for _ in range(len(df_n))], index=df_n.index)


# In[18]:


# Define a function to extract features
def extract_features(df, feature_types):
    for feature_type in feature_types:
        for i, row in df.iterrows():
            with open("tmp.fasta", "w") as f:
                f.write(">seq\n" + row["Sequence"])
            command = "/Users/anandnambigovindarajan/iFeature/iFeature.py --file tmp.fasta --type " + feature_type + " --out tmp.csv"
            subprocess.call(command, shell=True)
            if os.path.isfile("tmp.csv"):
                features_df = pd.read_csv("tmp.csv", sep='\t')
                features = [float(x) for x in features_df.values[0][1:]]
                if 'Features' in df.columns and df.at[i, "Features"]:
                    df.at[i, "Features"] += features
                else:
                    df.at[i, "Features"] = features
                os.remove("tmp.csv")
            os.remove("tmp.fasta")

# Define the feature types you want to extract
feature_types = ["AAC", "DPC", "CTDC", "CTDT"]


# In[19]:


# Extract features for the positive and negative data
extract_features(df, feature_types)
extract_features(df_n, feature_types)


# In[20]:


# Concatenate the positive and negative dataframes
df_all = pd.concat([df, df_n])


# In[23]:


# Split the data into a training set and a test set
X = df_all["Features"].tolist()
y = df_all["Label"].tolist()
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# In[44]:


from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression


# In[45]:


# Define the classifiers
classifiers = {
    'SVM': svm.SVC(),
    'KNN': KNeighborsClassifier(),
    'Decision Tree': DecisionTreeClassifier(),
    'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
    'Logistic Regression': LogisticRegression()
}


# In[46]:


# Train and evaluate each classifier
for name, clf in classifiers.items():
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    print(name)
    print(classification_report(y_test, y_pred))


# In[24]:


# Train a RandomForestClassifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)


# In[25]:


# Evaluate the model
y_pred = clf.predict(X_test)
print(classification_report(y_test, y_pred))


# In[ ]:


#######Prediction########


# In[32]:


# Assume that new_seq is your new sequence
new_sequence = "MDQEAGFMVNFINSYFIALGVLIGGALIGGLGAYLAGEPPLTAITKLANRLKIWALVAAIGGTFDAVYSFERGILEGNTRDIFKQLLLIISAMGGAQSGWLIISWLTQEHLSS"


# In[33]:


# Extract features from the sequence
def extract_features_from_sequence(sequence, feature_types):
    features = []
    with open("tmp.fasta", "w") as f:
        f.write(">seq\n" + sequence)
    for feature_type in feature_types:
        command = "/Users/anandnambigovindarajan/iFeature/iFeature.py --file tmp.fasta --type " + feature_type + " --out tmp.csv"
        subprocess.call(command, shell=True)
        if os.path.isfile("tmp.csv"):
            features_df = pd.read_csv("tmp.csv", sep='\t')
            features += [float(x) for x in features_df.values[0][1:]]
            os.remove("tmp.csv")
    os.remove("tmp.fasta")
    return features


# In[34]:


# Define the feature types you want to extract
feature_types = ["AAC", "DPC", "CTDC", "CTDT"]


# In[35]:


# Extract features
new_features = extract_features_from_sequence(new_sequence, feature_types)


# In[37]:


# Predict using the trained model
prediction = clf.predict([new_features])

if prediction[0] == 0:
    print("The sequence is predicted to be non-virulent")
else:
    print("The sequence is predicted to be virulent")


# In[38]:


# Assume that new_seq is your new sequence
new_sequence = "MNRREFLLNSTKTMFGTAALASFPLSIQKALAIDAKVESGTIQDVKHIVILTQENRSFDNYFGTLKGVRGFGDRFTIPMTEGRKVWEQYDANKKKVLPYHLDSRLGNAQRVTGTNHSWSDGQGAWDNGRMSDWVAHKQPQSMGYYKKQEVEYQFALANAFTICDAYHCAMHAGTNPNRKFIWTGTNGPTGAGVASVVNEFDGIGPSTEGYEWTTYPERLQQAGVTWKVYQNMPDNFTDNPLAGFKQYRRANEQSGQPVSNDTLICLAYDEKIDATQPLYKGIANTMPDGGFLGAFKADIAQGKLPQVSWLVAPATYSEHPGPSSPVQGAWYIQEVLNALTENTQVWSQTVLLVNFDENDGFFDHVPSPSAPSKDINGVVYGKTTLTDQQVSYEYFNHPAVATSKSQPETDGRVYGPGVRVPMYVISPWSRGGWVNSQVFDHTSILQFLEKRFGVQEPNISPYRRAVCGDLTTAFNFKTPNLLPVAELDGKKTKAEADAIRVAQELLPQVSVPSQQQFPQQEIGIRPSRALPYILHTSAKVDVTQKTVKLMFSNTGKQAAVFHVYNRLDLTAIPRRYMVEAGKQLDDAWNTINGQYDLWVLGPNGFHRAFKGNLSQANQTQALPEIRVCVEECDANLYLKVRHDGNKSVKLNVKANAYLPNKTWMIETNSSEKELVWDMSEFGGWYDFTVTLADDATFSRRFAGRIETQEDSISDPYMGYLES"


# In[39]:


# Extract features from the sequence
def extract_features_from_sequence(sequence, feature_types):
    features = []
    with open("tmp.fasta", "w") as f:
        f.write(">seq\n" + sequence)
    for feature_type in feature_types:
        command = "/Users/anandnambigovindarajan/iFeature/iFeature.py --file tmp.fasta --type " + feature_type + " --out tmp.csv"
        subprocess.call(command, shell=True)
        if os.path.isfile("tmp.csv"):
            features_df = pd.read_csv("tmp.csv", sep='\t')
            features += [float(x) for x in features_df.values[0][1:]]
            os.remove("tmp.csv")
    os.remove("tmp.fasta")
    return features


# In[40]:


# Define the feature types you want to extract
feature_types = ["AAC", "DPC", "CTDC", "CTDT"]


# In[41]:


# Extract features
new_features = extract_features_from_sequence(new_sequence, feature_types)


# In[42]:


# Predict using the trained model
prediction = clf.predict([new_features])

if prediction[0] == 0:
    print("The sequence is predicted to be non-virulent")
else:
    print("The sequence is predicted to be virulent")


# In[ ]:




