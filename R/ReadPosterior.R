
readPosterior <- function(files, par = "beta_snp", equal.par = F) {

  createCustomAwk <- function(is) {
    str <- paste('{print ', paste(paste("$", is, sep = ''),
                                  collapse = ','), "}", sep = '')
    return (str)
  }

  data <- vector(mode = "list",
                 length = length(files))

  for(f in 1:length(files)) {
    cat(f, ", ")

    if(file.exists(files[f]) == FALSE) {
      warning("File not found.")
      stop()
    }

    # mcmc.warmup
    awk.path.pipe <- pipe(paste("gawk -F , 'FNR == 10'",
                                files[f], sep = ' '))
    cols <- readLines(con = awk.path.pipe)
    close(con = awk.path.pipe)
    mcmc.warmup <- gsub(pattern = '\\#| |warmup|=',
                        replacement = '', x = cols)
    mcmc.warmup <- as.numeric(mcmc.warmup)


    # colnames
    awk.path.pipe <- pipe(paste("gawk -F , 'FNR == 26'",
                                files[f], sep = ' '))
    cols <- readLines(con = awk.path.pipe)
    close(con = awk.path.pipe)
    cols <- unlist(strsplit(x = cols, split = ','))

    if(equal.par == FALSE) {
      is <- which(regexpr(pattern = par, text = cols,
                          ignore.case = TRUE) != -1)
    }
    if(equal.par == TRUE) {
      is <- which(cols == par)
    }


    if(length(is) == 0) {
      stop("no such parameter found in sampling file.")
    }

    awk <- createCustomAwk(is)
    writeLines(text = awk, con = "awk_temp.awk")

    awk.path.pipe <- pipe(paste("gawk -F , -f awk_temp.awk",
                                files[f], "> temp.csv", sep = ' '))
    d <- readLines(con = awk.path.pipe)
    close(con = awk.path.pipe)

    awk.path.pipe <- pipe("gawk -F , '(NR > 25)' temp.csv > temp2.csv")
    d <- readLines(con = awk.path.pipe)
    close(con = awk.path.pipe)

    # here take only samples (-warmup)
    d <- read.table(file = "temp2.csv", header = TRUE,
                    as.is = TRUE, sep = ' ')
    d <- d[complete.cases(d), ]

    if(is.vector(d)) {
      d <- matrix(data = d, ncol = 1)
    }

    d <- d[-(1:mcmc.warmup), ]
    data[[f]] <- d

    # do some cleanup before next file
    file.remove(c("awk_temp.awk", "temp.csv", "temp2.csv"))
  }

  return (data)
}
