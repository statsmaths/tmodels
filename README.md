# tmodels: Apply Statistical Tests with Consistent Formula Interface

Check contingency table tests

```{r}
library(tmodels)

df <- data.frame(x = sample(c("H", "T"), 100, replace=TRUE),
                 y = sample(c("A", "B"), 100, replace=TRUE))
df$y[df$x == "H"] <- sample(c("A", "A", "A", "B"), sum(df$x == "H"), replace=TRUE)

tmod_z_test_prop(y ~ x, data = df)
tmod_odds_ratio(y ~ x, data = df)
tmod_fisher_test(y ~ x, data = df)
tmod_chi_squared_test(y ~ x, data = df)

df <- data.frame(x = sample(c("H", "T", "G"), 100, replace=TRUE),
                 y = sample(c("A", "B"), 100, replace=TRUE))

tmod_odds_ratio(y ~ x, data = df)
tmod_fisher_test(y ~ x, data = df)
tmod_chi_squared_test(y ~ x, data = df)

df <- data.frame(x = sample(c("H", "T", "G"), 100, replace=TRUE),
                 y = sample(c("A", "B", "C"), 100, replace=TRUE))

tmod_odds_ratio(y ~ x, data = df)
tmod_fisher_test(y ~ x, data = df)
tmod_chi_squared_test(y ~ x, data = df)
```

check correlation tests

```{r}
library(tmodels)

df <- data.frame(x = runif(100),
                 y = runif(100))
df$y <- df$x + df$y * 0.1

tmod_pearson_correlation_test(y ~ x, data = df)
tmod_spearman_correlation_test(y ~ x, data = df)
tmod_kendall_correlation_test(y ~ x, data = df)
```

check difference in means tests (2 groups)

```{r}
library(tmodels)

df <- data.frame(x = sample(c("H", "T"), 100, replace=TRUE),
                 y = runif(100))
df$y[df$x == "H"] <- df$y[df$x == "H"] + 0.2

t.test(y ~ x, data = df)
tmod_t_test(y ~ x, data = df)

wilcox.test(y ~ x, data = df)
tmod_mann_whitney_test(y ~ x, data = df)
```

check difference in means tests (3+ groups)

```{r}
library(tmodels)

df <- data.frame(x = sample(c("H", "T", "C"), 100, replace=TRUE),
                 y = runif(100))
df$y[df$x == "H"] <- df$y[df$x == "H"] + 0.2

anova(lm(y ~ x, data = df))
tmod_one_way_anova_test(y ~ x, data = df)

kruskal.test(y ~ x, data = df)
tmod_kruskal_wallis_test(y ~ x, data = df)
```

