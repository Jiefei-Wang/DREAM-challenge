library(synapser)
synLogin()
mySubmission <- File("/path/to/submission.csv", parentId="syn12345")
mySub <- synStore(mySubmission)
#The evaluationId HAS to be a string here or there will be an error
submission <- synSubmit(evaluation, mySub, name="Our Final Answer", team="Blue Team") 