library(arrow)
###############################################################################
df_a <- data.frame(id = 1:5, value = rnorm(5), subset = "a")
df_b <- data.frame(id = 6:10, value = rnorm(5), subset = "b")
df_c <- data.frame(id = 11:15, value = rnorm(5), subset = "c")
ds_dir <- "mini-dataset"
dir.create(ds_dir)
ds_dir_a <- file.path(ds_dir, "subset=a")
ds_dir_b <- file.path(ds_dir, "subset=b")
ds_dir_c <- file.path(ds_dir, "subset=c")

dir.create(ds_dir_a)
dir.create(ds_dir_b)
dir.create(ds_dir_c)
write_parquet(df_a, file.path(ds_dir_a, "part-0.parquet"))
write_parquet(df_b, file.path(ds_dir_b, "part-0.parquet"))
write_parquet(df_c, file.path(ds_dir_c, "part-0.parquet"))
list.files(ds_dir, recursive = TRUE)
ds <- open_dataset(ds_dir)
glimpse(ds)
ds$files
ds$schema
scan <- Scanner$create(dataset = ds)
scan$ToTable()
ds.inmemory <- as.data.frame(scan$ToTable())
ds.batch <- lapply(scan$ScanBatches(), as.data.frame)

###############################################################################
ds <- open_dataset("nyc-taxi")
ds <- open_dataset("nyc-taxi", partitioning = c("year", "month"))
###############################################################################
####### logistic-simu-dataset
library(arrow)
simu_mzNormal <- function(N, beta0, corr, id_range){
  d <- length(beta0) - 1
  sigmax  <- matrix(corr, d, d) + diag(1-corr, d)
  X <- mvtnorm::rmvnorm(N, rep(0, d), sigmax)
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- rbinom(N, 1, P)

  return(data.frame(
    id = id_range,
    Y = Y,
    X = X,
    subset = as.integer(id_range[1])
  ))
}
ds_dir <- "logistic-simu-dataset"
if (!dir.exists(ds_dir)) dir.create(ds_dir)
N <- 3e7
# Loop
for (i in 1:10) {
  df <-
    simu_mzNormal(
      N = N,
      beta0 = rep(0.5, 7),
      corr = 0.5,
      id_range = c((1:N) + (i - 1) * N)
    )

  subset_dir <- file.path(ds_dir, paste0("subset=", i))
  if (!dir.exists(subset_dir)) dir.create(subset_dir)
  write_parquet(df, file.path(subset_dir, paste0("part-0.parquet")))
  rm(df)
  gc()
  # gc(verbose = TRUE)
}
ds <- open_dataset(ds_dir)
model <- lm(Y ~ X.1, data = ds)
summary(model)


# arrow searches the dataset folder looking for appropriate files,
# but does not load the contents of those files. Paths to these files are stored
# in an active binding ds$files:
ds$files
ds$schema

scan <- Scanner$create(dataset = ds)
# Calling the ToTable() method will materialize the Dataset (on-disk) as a Table (in-memory):
ds.scan <- scan$ToTable()
ds.in.memory <- as.data.frame(ds.scan)
object.size(ds.scan)
object.size(ds)
object.size(ds.in.memory)
model <- lm(Y ~ X.1, data = ds.scan)
model <- lm(Y ~ X.1, data = ds.in.memory)
