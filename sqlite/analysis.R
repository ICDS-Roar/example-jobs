install.packages("DBI")
install.packages("RSQLite")
library(DBI)
library(tidyverse)

# Create a new database (or open an existing one, if it already exists)
con <- dbConnect(RSQLite::SQLite(), "state.db")

# Add general state information - from `datasets` package
data(state)
state_data <- tibble(name=state.name, abbr=state.abb, region=state.region)
dbWriteTable(con, "states", state_data)

# Add some income and rental data - from `tidyr` package
data(us_rent_income)
names(us_rent_income) <- c("ID", "name", "variable", "estimate", "error")

us_monthly_rents <- us_rent_income %>% 
  filter(variable == "rent") %>% 
  select(name, estimate, error) %>%
  rename(rental_err = error, ave_mo_rent = estimate)

us_yearly_income <- us_rent_income %>% 
  filter(variable == "income") %>% 
  select(name, estimate, error) %>%
  rename(income_err = error, ave_yr_income = estimate)

us_rent_income <- left_join(us_monthly_rents, us_yearly_income)
dbWriteTable(con, "rental_costs", us_rent_income)

# Query the database using SQL syntax
query <- "SELECT s.region, s.name, ave_mo_rent/(ave_yr_income/12) as 'rent_income'
          FROM (SELECT name, region FROM states) s
          JOIN rental_costs c on s.name = c.name
          ORDER BY rent_income DESC;"
results <- dbGetQuery(con, query)


# create plot from data
ggplot(results, aes(fill=region, y=rent_income, x=reorder(name, rent_income))) +
  geom_bar(stat="identity", position="dodge") + 
  labs(title = "Average Rent to Income Ratio", x="State", fill="Region") + 
  scale_x_discrete(guide=guide_axis(angle=70))

