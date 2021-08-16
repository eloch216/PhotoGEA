# Try to apply as.numeric, but if any warnings or errors occur, just use the
# original column values
try_as_numeric <- function(x) {
    tryCatch(
        {
            # Code to be executed initially
            as.numeric(x)
        },
        error=function(cond) {
            # Code for handling errors
            x
        },
        warning=function(cond) {
            # Code for handling warnings
            x
        }
    )
}
