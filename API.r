library(plumber)
library(attempt)

api = plumber::plumb("/home/daqop/Desktop/cours/stage/R/repo.r")
api$run(host="127.0.0.1", port=8048)