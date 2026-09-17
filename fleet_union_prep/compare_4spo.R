# Field-by-field serialized comparison of a rebuilt camera cache against the
# cached 4SPO of record. grid100 (dead terra pointer) and camera_prep (new
# provenance field) are excluded. Exit status 1 on any difference.
args <- commandArgs(trailingOnly = TRUE)
A <- readRDS(args[1]); B <- readRDS(args[2])
flds <- setdiff(union(names(A), names(B)), c("grid100", "camera_prep"))
bad <- 0
ser <- function(x) serialize(x, NULL, version = 3)
for (f in flds) {
  if (is.list(A[[f]]) && !is.data.frame(A[[f]])) {
    for (g in union(names(A[[f]]), names(B[[f]]))) {
      same <- identical(ser(A[[f]][[g]]), ser(B[[f]][[g]]))
      if (!same) { bad <- bad + 1; cat("   DIFF", f, "$", g, "\n") }
    }
  } else {
    same <- identical(ser(A[[f]]), ser(B[[f]]))
    if (!same) {
      bad <- bad + 1
      num <- tryCatch(all.equal(as.data.frame(A[[f]]), as.data.frame(B[[f]]), check.attributes = FALSE), error = function(e) conditionMessage(e))
      cat("   DIFF", f, ":", paste(head(num, 3), collapse = " | "), "\n")
    }
  }
}
cat(if (bad == 0) "   CAMERA CACHE IDENTICAL\n" else paste("   CAMERA CACHE DIFFERS in", bad, "field(s)\n"))
quit(status = if (bad == 0) 0 else 1)
