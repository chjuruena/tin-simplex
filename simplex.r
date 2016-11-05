E1 <- function (x1, x2) 7 * x1 + 11 * x2 + -77;
E2 <- function (x1, x2) 10 * x1 + 8 * x2 + -80;
E3 <- function (x1, x2) 1 * x1 + 0 * x2 + -9;
E4 <- function (x1, x2) 0 * x1 + 1 * x2 + -6;
E5 <- function (x1, x2) -150 * x1 + -175 * x2 + 0;
f <- list(E1, E2, E3, E4, E5);


# a = c(7,11,1,0,0,0,0,77);
# b = c(10,8,0,1,0,0,0,80);
# c = c(1,0,0,0,1,0,0,9);
# d = c(0,1,0,0,0,1,0,6);
# e = c(-150, -175, 0,0,0,0,1,0);

# mat = rbind(a,b,c,d,e);

Simplex <- function(f, verbose){
  
  mat = AugCoeffMatrix(f)$augcoeffmatrix;
  var = AugCoeffMatrix(f)$variables;
  row = nrow(mat);
  mat = cbind (mat[,1], mat[,2], diag(x = 1, row ), mat[,3]);
  
  col = ncol(mat);
  
  
  
  slack = c(paste("s", 1:row)); #setting slack variables after setting the new row and col ofthe matrix
  
  # label = c("x1","x2", "S1", "S2", "S3", "S4", "Z", "S");
  label = c(var, slack, "S");
  dimnames(mat)=list(c(1:(row)), label);
  
  print(mat)
  iter = 1;
  solutionSet = c();
  # stop only if all the last row is positive already
  while( min(mat[row, ]) < 0){
    
    solutionSet = c();
    #get max value on the last row
    pivotcol = which.min((mat[row, ]))
    
    # get row of min solution/x
    
    temp = c();
    for (i in 1:row){			
      x = as.numeric(mat[i, col]/mat[i,pivotcol]);
      if(as.numeric(mat[i, col] == 0) || ( x == Inf) || ( x < 0)  ){
        temp = c(temp, Inf);
      }
      else {
        temp = c(temp, x);
      }
      
    }
    pivotrow = as.numeric(which.min(abs(temp)));
    # test ratio - smallest positive
    
    
    
    #Normalization
    mat[pivotrow, ] = round(mat[pivotrow, ] / as.numeric(mat[pivotrow, pivotcol]), digits = 4); #Normalization 
    
    # set dim names for the matrix
    # static only
    # label = c("x1","x2", "S1", "S2", "S3", "S4", "Z", "S");
    # 		dimnames(mat)=list(c(1:(row)), label);
    
    for(j in 1:row){
      #multiplier = mat[j,i]/mat[i,i];
      
      if(as.numeric(mat[j,pivotcol]) == 0 || (j == pivotrow)) {
        next;
      }
      
      
      normalized =  round(mat[pivotrow,] * as.numeric(mat[j,pivotcol]), digits = 4); #current row * term sa evalrow
      mat[j,] = round(mat[j,] - normalized, digits = 4); #evalrow=evalrow-normalized   
      
      if(verbose){
        print(paste("Resulting Matrix after new eval_row", j, "iteration", iter),quote = FALSE);
        print(mat);
      }
      
    }
    
    
    solutionSet = c();
    
    # solution set
    for(j in 1:col){
      cnt = 0;
      for(k in 1:row){
        if (mat[k,j] == 0){
          cnt = cnt + 1;
        }
        if(mat[k,j] == 1){
          temp = mat[k,col];
        }
      }
      if ( (max(mat[ ,j]) == 1 ) && (cnt<=row || (cnt==(row-1)))){
        solutionSet = c(solutionSet, temp);
      }else{
        solutionSet = c(solutionSet, 0);
      }
    }
    names(solutionSet)=label;
    
    
    
    if(verbose){
      
      print(paste("******************************************"),quote = FALSE);	
      print(paste("Resulting Matrix after Iteration", iter),quote = FALSE);
      print(mat);
      print(paste("Solution Set: "),quote = FALSE);	
      print(solutionSet)
      print(paste("******************************************"),quote = FALSE);	
      
    }
    
    iter=iter+1;
    if(min((mat[row, ]))>0) break; #if the smallest in the last row is  already positive then stop
  }
  if(verbose){
    print(paste("Resulting Matrix for Simplex"),quote = FALSE);
    print(mat);
    print(paste("Solution Set: "),quote = FALSE);	
    print(solutionSet)
  }
  
  
  result = list(mat = mat, solutionSet = solutionSet, x1 = solutionSet[1], x2 = solutionSet[2]);
  
  return (result);
  
}

AugCoeffMatrix <- function(mylist){
  i = 1;
  l = c();
  #deparsing the list per item
  len = length(mylist);
  for (i in 1:length(mylist)){
    temp = deparse((mylist[[i]]));
    l = c(l , temp);
  }
  
  i=1;
  #from the deparsed list, traverse list and split  
  for (i in 1:length(l)){
    if(i %% 2 == 0){
      temp = l[i];
      #split string by space
      temp2 = unlist(strsplit(temp, "[[:space:]]"));
      j=1;
      #using grep to find the values by regex
      vars = grep("x+", temp2, perl=TRUE, value=TRUE);
      line = grep("^[-+]?[0-9]*\\.?[0-9]+", temp2, perl=TRUE, value=TRUE);
      
      newval = c();
      #order the variables
      for (x in 1:length(vars)){
        index = substring(vars[x],2)
        index = as.numeric(index)
        
        if (x == index ){
          newval = c(newval, line[x])
        }
        else{
          #newval[index+1] = line[x]
          newval[index] = line[x]
          
        }
      }
      newval[length(vars)+1] = line[length(vars)+1];
      
      
      #converts to a number
      newval = as.numeric(newval)
      linelen = length(newval);
      newval[linelen] = newval[linelen]*-1;
      
      if (i == 2){
        mat = newval;
      }else{
        mat = rbind(mat, newval);
      }
    }
  }
  #set dimension names
  dimnames(mat)=list(c(1:len), c(vars, "RHS"));
  result = list (variables=c(vars),augcoeffmatrix=mat);
  return (result);
}

