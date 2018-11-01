# Lecture 9: Unsupervised Learning Analysis of Cancer Cells

# Import the data 
url <- "https://bioboot.github.io/bimm143_S18/class-material/WisconsinCancer.csv"
wisc.df <- read.csv(url)
# Use as.matrix() to convert the other features (i.e. columns) of the data 
# (in columns 3 through 32) to a matrix. Store this in a variable called wisc.data.
wisc.data <- as.matrix(wisc.df[3:32])
#Assign the row names of wisc.data the values currently 
#contained in the id column of wisc.df.
row.names(wisc.data) <- wisc.df$id
#Finally, setup a separate new vector called diagnosis to be 1 if 
#a diagnosis is malignant ("M") and 0 otherwise.
diagnosis <- as.numeric (wisc.df$diagnosis == "M")
