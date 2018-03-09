#Specify the state transition function:
#WARNING: always use arguments x and k when specifying the GGfunction
state_transition <- function (x, k){
    x + (A%*%x)*dt
}