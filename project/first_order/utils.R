denom = function(x,b,r,t,k){
  return(r + b*abs(x = t)**k)
}

loss_fn = function(x,y,a,b,r,k,t){
  return(sum(y - a*(x-t)/denom(x,b,r,t,k)**(1/k))**2)
}

grad_a = function(x,y,a,b,r,k,t){
  return(sum(-2*(x-t)/denom(x,b,r,t,k)**(1/k)*(y - a*(x-t)/denom(x,b,r,t,k)**(1/k))))
}

grad_r = function(x,y,a,b,r,k,t){
  first_chunk = 2*y*a*(x-t)/k*denom(x,b,r,t,k)*(-1/k-1)
  second_chunk = -2*a*(x-t)**2/k*denom(x,b,r,k,t)**(-2/k - 1)
  return(sum(first_chunk + second_chunk))
}

grad_b = function(x,y,a,b,r,k,t){
  first_chunk = 2*a*(x-t)*abs(x-t)**k*denom(a,b,r,t,k)**(-1/k -1)
  second_chunk = y - a*(x-t)*denom(a,b,r,t,k)**(-1/k)
  return(sum(first_chunk*second_chunk)/k)
}

grad_k = function(x,y,a,b,r,k,t){
  first_chunk = y - a*(x-t)/denom(x,b,r,t,k)**(1/k)
  second_chunk = -a*(x-t)*denom(x,b,r,t,k)**(-1/k)*log(denom(x,b,r,t,k))/k**2
  third_chunk = (denom(x,b,r,t,k) - r)*log(abs(x - t))/(k*denom(x,b,r,t,k))
  return(sum(2*first_chunk*(second_chunk - third_chunk)))
}

grad = function(x,y,a,b,r,k,t){
  return(c(grad_a(x,y,a,b,r,k,t),
           grad_b(x,y,a,b,r,k,t),
           grad_r(x,y,a,b,r,k,t),
           grad_k(x,y,a,b,r,k,t)))
}

update_x = function(x,y,a,b,r,k,t,gamma,delta){
  x_next = c(a,b,r,k) - gamma*grad(x,y,a,b,r,k,t) + alpha*delta
  return(x_next)
}

