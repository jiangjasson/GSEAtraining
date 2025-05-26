
#' A REPL in a local environment
#' 
#' @param parent Parent environment.
#' 
#' @export
interactive_local = function(parent = parent.frame()) {
	env = new.env(parent = parent)

	read_code = function() {

		exit_local = function() {

		}
		buffer = ""
		repeat{
			cat(if (nzchar(buffer)) "local + " else "local > ")  # REPL prompt
			line = readLines(stdin(), n = 1)
			buffer = paste(buffer, line, sep = "\n")

			is_complete = tryCatch({
				parse(text = buffer)
				TRUE
			}, error = function(e) {
				if (grepl("unexpected end of input", e$message)) {
					FALSE  # Incomplete code
				} else {
					message("Parse error: ", e$message)
					buffer <<- ""  # Reset on error
					FALSE
				}
			})

			if(is_complete) {
				break
			}
		}
		buffer
	}

	repeat {
		code = read_code()

		if(grepl("exit_local()", code)) {
			break
		}

		exprs = parse(text = code)
		for (expr in exprs) {
			result = try(eval(expr, envir = env), silent = TRUE)
				if (inherits(result, "try-error")) {
				cat("Error: ", result, "\n")
			} else if (!is.null(result)) {
				print(result)
			}
		}
	}

}
