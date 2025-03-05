sum_list = function(list1, list2){
  list_sum = list()
  for (key in names(list1)) {
    list_sum[[key]] <- list1[[key]] + list2[[key]]
  }
  return(list_sum)
}

sum_div = function(list, constant){
  list_div = list()
  for (key in names(list)) {
    list_div[[key]] <- list[[key]] / constant
  }
  return(list_div)
}