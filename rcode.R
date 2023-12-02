state <- c("Abia.1",        "Abuja.1" , "Adamawa.1",     "Akwa Ibom.1" ,  "Anambra.1"  ,   "Bauchi.1",      "Bayelsa.1",    
           "Benue.1",       "Borno.1" ,      "Cross River.1", "Delta.1",       "Ebonyi.1",      "Edo.1",         "Ekiti.1",      
           "Enugu.1"   ,    "Gombe.1"  ,     "Imo.1",         "Jigawa.1",      "Kaduna.1",    "Kano.1" ,       "Katsina.1",    
           "Kebbi.1",       "Kogi.1"  ,      "Kwara.1"  ,     "Lagos.1" ,    "Nassarawa.1" , 
           "Niger.1",       "Ogun.1" ,       "Ondo.1" ,       "Osun.1" ,       "Oyo.1"    ,     "Plateau.1"   ,  "Rivers.1",     
           "Sokoto.1",      "Taraba.1",      "Yobe.1",        "Zamfara.1")
numBUGS= c(22,37,8,1,2,3,31,5,6,7,23,32,4,33,24,34,9,25,10,11,12,26,27,13,14,35,15,16,17,28,18,19,20,21,29,30,36)

nstates = 2 #length(state)

W = matrix(0, nrow = nstates, ncol = nstates)

neighbour_list = list(c(22, 20, 7, 1, 9, 2, 24, 32), # needs to be at the same order as numBUGS 
                      c(37, 15, 10, 35, 27)) 

for(i in 1:nstates){
  where = NULL
  num_neighbours = length(neighbour_list[[i]][-1])
  
  for(j in 1:nstates){
    for(k in 1:num_neighbours){
      if(numBUGS[j] == neighbour_list[[i]][k + 1]){
        where[k] = j
      }
    }
  }
  W[i,where] = 1
}

D = diag(rowSums(W))