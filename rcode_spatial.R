state <- c("Abia.1",        "Abuja.1" , "Adamawa.1",     "Akwa Ibom.1" ,  "Anambra.1"  ,   "Bauchi.1",      "Bayelsa.1",    
           "Benue.1",       "Borno.1" ,      "Cross River.1", "Delta.1",       "Ebonyi.1",      "Edo.1",         "Ekiti.1",      
           "Enugu.1"   ,    "Gombe.1"  ,     "Imo.1",         "Jigawa.1",      "Kaduna.1",    "Kano.1" ,       "Katsina.1",    
           "Kebbi.1",       "Kogi.1"  ,      "Kwara.1"  ,     "Lagos.1" ,    "Nassarawa.1" , 
           "Niger.1",       "Ogun.1" ,       "Ondo.1" ,       "Osun.1" ,       "Oyo.1"    ,     "Plateau.1"   ,  "Rivers.1",     
           "Sokoto.1",      "Taraba.1",      "Yobe.1",        "Zamfara.1")
numBUGS= c(22,37,8,1,2,3,31,5,6,7,23,32,4,33,24,34,9,25,10,11,12,26,27,13,14,35,15,16,17,28,18,19,20,21,29,30,36)

nstates = 37 #length(state)

W = matrix(0, nrow = nstates, ncol = nstates)

neighbour_list = list(c(22, 20, 7, 1, 9, 2, 24, 32), # needs to be at the same order as numBUGS 
                      c(37, 15, 10, 35, 27), c(8, 29, 34, 6), c(1, 20, 22, 7), c(2, 27, 24, 22, 9, 23, 4, 20), 
                      c(3, 25, 11, 19, 29, 34, 10, 30), c(31, 20, 23), c(5, 35, 29, 27, 24, 32, 7), c(6, 30, 34, 8), 
                      c(7, 5, 32, 22, 1), c(23, 17, 4, 20, 2, 31), c(32, 24, 5, 7, 22), c(4, 17, 23, 27, 2), 
                      c(33, 28, 27, 13, 27), c(24, 27, 5, 32, 2, 22), c(34, 3, 30, 6, 8, 29), c(9, 20, 22, 2), 
                      c(25, 11, 3, 30, 12), c(10, 15, 36, 12, 11, 3, 19, 35, 37), c(11, 25, 12, 3, 10), 
                      c(12, 36, 10, 11, 25), c(26, 15, 36, 21), c(27, 13, 15, 37, 35, 5, 24, 2, 33, 17, 4), 
                      c(13, 15, 18, 28, 33, 27), c(14, 16), c(35, 10, 37, 19, 29, 5, 11), c(15, 13, 27, 37, 10, 36, 26), 
                      c(16, 14, 28, 18, 17), c(17, 16, 28, 33, 4, 23, 27), c(28, 18, 16, 17, 33, 13), c(18, 13, 28, 16), 
                      c(19, 10, 3, 35, 29), c(20, 23, 31, 9, 22, 2, 1), c(21, 26, 36), c(29, 5, 35, 19, 3, 34, 8), 
                      c(30, 3, 34, 6), c(36, 26, 15, 21, 12, 10)) 

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

D_W_inv = solve(D-W) # enter on stan as data

##### Shapefile####
state <- c("Abia.1",        "Abuja.1" , "Adamawa.1",     "Akwa Ibom.1" ,  "Anambra.1"  ,   "Bauchi.1",      "Bayelsa.1",    
           "Benue.1",       "Borno.1" ,      "Cross River.1", "Delta.1",       "Ebonyi.1",      "Edo.1",         "Ekiti.1",      
           "Enugu.1"   ,    "Gombe.1"  ,     "Imo.1",         "Jigawa.1",      "Kaduna.1",    "Kano.1" ,       "Katsina.1",    
           "Kebbi.1",       "Kogi.1"  ,      "Kwara.1"  ,     "Lagos.1" ,    "Nassarawa.1" , 
           "Niger.1",       "Ogun.1" ,       "Ondo.1" ,       "Osun.1" ,       "Oyo.1"    ,     "Plateau.1"   ,  "Rivers.1",     
           "Sokoto.1",      "Taraba.1",      "Yobe.1",        "Zamfara.1")
numBUGS= c(22,37,8,1,2,3,31,5,6,7,23,32,4,33,24,34,9,25,10,11,12,26,27,13,14,35,15,16,17,28,18,19,20,21,29,30,36)
data.frame(state,numBUGS)
adjBUGS <- data.frame(region=state,
                      est=secd$sstate[1:37][numBUGS])
#for (i in 1:37) {
#  shp.map <-s(node1[i,1],node2[i,2])
# node1 <- 
#}
#s <- function(region,est){
# shp.map$estimate[which(shp.map$group==as.character(region))] <- est
#  return(shp.map)
#}
mungeCARdata4stan = function(adjBUGS,numBUGS) {
  k = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:k) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N"=k,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}