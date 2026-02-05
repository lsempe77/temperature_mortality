library(microdatasus)

res <- fetch_datasus(year_start = 2021, 
                     year_end = 2021,
                     timeout= 500,
                     information_system = "SIM-DO")

res_proc <-process_sim(res)

write.csv(res_proc,"DO21OPEN.csv")


