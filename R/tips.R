## how to add maxk k in sm?
maxk=20
mgcv::s

mgcv::gam

getAnywhere("interpret.gam0")

formular = numvisits ~ children + race + maritalstat + 
  s(age, bs = "ps", k = maxk)+ 
  s(income1000, bs = "ps", k = maxk)+
  s(access, bs = "ps", k = maxk)+
  s(pc1times1000, bs = "ps", k = maxk)
mgcv::interpret.gam(formular)
tf= terms.formula(formular, specials = "s")
tf
attr(tf, "specials")
attr(tf, "term.labels")
attr(tf, "term.labels")[attr(tf, "specials")[[1]] - 1][1]
parse(text = paste("mgcv::", 
                   "s(age, bs = \"ps\", k = maxk)", sep = ""))
x = try(eval(parse(text = paste("mgcv::", 
                            "s(age, bs = \"ps\", k = maxk)", sep = ""))))
x


function (..., k = -1, fx = FALSE, bs = "tp", m = NA, by = NA, 
          xt = NULL, id = NULL, sp = NULL, pc = NULL) 
{
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  if (by.var == ".") 
    stop("by=. not allowed")
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (term[1] == ".") 
    stop("s(.) not supported.")
  if (d > 1) 
    for (i in 2:d) {
      term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
      if (term[i] == ".") 
        stop("s(.) not yet supported.")
    }
  for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])), 
                                 "term.labels")
  k.new <- round(k)
  if (all.equal(k.new, k) != TRUE) {
    warning("argument k of s() should be integer and has been rounded")
  }
  k <- k.new
  if (length(unique(term)) != d) 
    stop("Repeated variables as arguments of a smooth are not permitted")
  full.call <- paste("s(", term[1], sep = "")
  if (d > 1) 
    for (i in 2:d) full.call <- paste(full.call, ",", 
                                      term[i], sep = "")
  label <- paste(full.call, ")", sep = "")
  if (!is.null(id)) {
    if (length(id) > 1) {
      id <- id[1]
      warning("only first element of `id' used")
    }
    id <- as.character(id)
  }
  ret <- list(term = term, bs.dim = k, fixed = fx, dim = d, 
              p.order = m, by = by.var, label = label, xt = xt, id = id, 
              sp = sp)
  if (!is.null(pc)) {
    if (length(pc) < d) 
      stop("supply a value for each variable for a point constraint")
    if (!is.list(pc)) 
      pc <- as.list(pc)
    if (is.null(names(pc))) 
      names(pc) <- unlist(lapply(vars, all.vars))
    ret$point.con <- pc
  }
  class(ret) <- paste(bs, ".smooth.spec", sep = "")
  ret
}