library(foreach)
library(doParallel)

if(exists("cl")){
  try(stopCluster(cl=cl),silent=T)
}

master="192.168.1.48"

machineAddresses <- list(
  list(host='localhost',ncore=10),
  list(host='192.168.1.230',user='jeff',rscript="Rscript",rshcmd="plink -pw qwer",
        ncore=40)
  # ,list(host='192.168.1.220',user='jeff',rscript="Rscript",rshcmd="plink -pw qwer",
  #      ncore=4)
  # ,list(host='192.168.1.217',user='jeff',rscript="Rscript",rshcmd="plink -pw qwer",
  #      ncore=4)
  # ,list(host='192.168.1.214',user='jeff',rscript="Rscript",rshcmd="plink -pw qwer",
  # ncore=8)
)

# machineAddresses <- list(
# list(host='192.168.0.102',user='jeff',rscript="Rscript",rshcmd="plink -pw qwer",
#         ncore=1)
# )



spec <- lapply(machineAddresses,
               function(machine) {
                 if(machine$host=="localhost")
                   rep(list(list(host=machine$host)),
                       machine$ncore)
                 else
                   rep(list(list(master=master,
                                 host=machine$host,
                                 user=machine$user,
                                 rscript=machine$rscript,
                                 rshcmd=machine$rshcmd
                                 )),
                       machine$ncore)
               })
spec <- unlist(spec,recursive=FALSE)



cl <- makePSOCKcluster(spec,manual = F)
# cl=makeCluster(11)
doParallel::registerDoParallel(cl)