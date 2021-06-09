condition <- function(subclass, message, call = sys.call(-1), ...) {
  structure(
    class = c(subclass, "condition"),
    list(message = message, call = call),
    ...
  )
}

userErrorCls <- "SVDFunctionsException"
userError <- function(msg) {
  e <- condition(c(userErrorCls), msg)
  stop(e)
}

