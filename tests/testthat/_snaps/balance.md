# dropping an unused exposure level announces itself

    Code
      fit <- balance(data, exposure, c(x1, x2), method = bw_entropy())
    Message
      i Treating `.exposure` as binary
      i Dropping unused exposure level "2".

