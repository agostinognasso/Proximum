# The loans dataset -----------------------------------------------------------
#
# `vignettes/credit-scoring-case.Rmd` asks three questions of a credit
# portfolio, and a synthetic portfolio is only worth shipping if all three have
# an answer that is in the data rather than in the telling. Real lending data
# that can be redistributed does not exist: the public credit panels are either
# licensed or stripped of the protected attributes the third question is about.
# So the portfolio is generated, and generated to a specification:
#
#   1. Risk segments the forest can find. Default depends on the score, the
#      utilisation and the debt-to-income ratio, and on an interaction between
#      the score and the utilisation, so that the proximity has blocks in it
#      rather than a single gradient.
#   2. A shift between vintages. Utilisation carries three times the weight in
#      the 2022 vintage that it carries in 2019, so "does the forest reorganise
#      its view of the portfolio from one year to the next" is a question with
#      a real answer rather than sampling noise.
#   3. A protected attribute the model never sees. `applicant_group` shifts the
#      income and the score, and enters the default probability nowhere at all.
#      A forest fitted without it will still separate the groups, because the
#      predictors carry them, and `permanova()` on the proximity will say how
#      much. That is the fairness diagnostic the vignette is for, and it needs
#      an attribute that is genuinely absent from the model and genuinely
#      present in the representation.
#
#   Rscript inst/data-raw/loans.R
#
# Writes data/loans.rda. Takes a second.

set.seed(20260908)

n <- 2000L
vintage <- factor(sample(2019:2022, n, replace = TRUE, prob = c(.22, .26, .26, .26)))

# The protected attribute, and the two things it moves. The gap is deliberate
# and is the point of question 3: it is the size of gap that a lender would
# call a difference in the applicant pool rather than a difference in the
# model.
applicant_group <- factor(sample(c("A", "B"), n, replace = TRUE, prob = c(.72, .28)))
group_shift <- ifelse(applicant_group == "B", 1L, 0L)

income <- round(exp(rnorm(n, mean = 10.9 - 0.22 * group_shift, sd = 0.45)))
score <- round(pmin(850, pmax(300,
  rnorm(n, mean = 690 - 28 * group_shift, sd = 62))))

employment <- pmin(40, round(rexp(n, rate = 1 / 6)))
home <- factor(sample(c("rent", "mortgage", "own"), n, replace = TRUE,
                      prob = c(.42, .45, .13)),
               levels = c("rent", "mortgage", "own"))
purpose <- factor(sample(
  c("debt consolidation", "home improvement", "major purchase", "medical",
    "small business"),
  n, replace = TRUE, prob = c(.52, .16, .14, .09, .09)
))
term <- factor(sample(c("36 months", "60 months"), n, replace = TRUE,
                      prob = c(.71, .29)))

amount <- round(pmin(40000, pmax(1000,
  exp(rnorm(n, mean = 9.3 + 0.35 * log(income / 30000), sd = 0.5)))) / 25) * 25
dti <- round(pmin(0.6, pmax(0.01,
  rnorm(n, mean = 0.19 + 0.9 * amount / income, sd = 0.06))), 3)
utilisation <- round(pmin(1, pmax(0,
  rnorm(n, mean = 0.55 - 0.0009 * (score - 690), sd = 0.22))), 3)
delinquencies <- rpois(n, lambda = pmax(0.02, 0.9 - 0.0035 * (score - 560)))

# Question 1 lives in the interaction and question 2 in the vintage weight. A
# thin borrower with a high utilisation is the segment; in 2019 the utilisation
# barely matters and by 2022 it dominates.
vintage_weight <- c("2019" = 0.6, "2020" = 1.0, "2021" = 1.5, "2022" = 1.9)
risk <- -4.00 -
  0.011 * (score - 690) +
  4.2 * (dti - 0.22) +
  vintage_weight[as.character(vintage)] * 1.5 * (utilisation - 0.55) +
  0.34 * delinquencies +
  0.55 * (term == "60 months") +
  0.22 * (home == "rent") -
  0.020 * employment +
  2.6 * (score < 660) * (utilisation > 0.7)

default <- factor(ifelse(rbinom(n, 1L, stats::plogis(risk)) == 1L, "yes", "no"),
                  levels = c("no", "yes"))

loans <- data.frame(
  default, amount, income, dti, term, employment, score, utilisation,
  delinquencies, home, purpose, vintage, applicant_group,
  stringsAsFactors = FALSE
)
loans <- loans[order(loans$vintage), ]
rownames(loans) <- NULL

# --- what the specification actually produced --------------------------------

cat("rows:", nrow(loans), " default rate:",
    sprintf("%.1f%%", 100 * mean(loans$default == "yes")), "\n")
cat("default rate by vintage:\n")
print(round(tapply(loans$default == "yes", loans$vintage, mean), 3))
cat("default rate by applicant group:\n")
print(round(tapply(loans$default == "yes", loans$applicant_group, mean), 3))
cat("the group is absent from the model and present in the predictors:\n")
print(round(c(
  income = diff(tapply(loans$income, loans$applicant_group, median)),
  score = diff(tapply(loans$score, loans$applicant_group, median))
), 1))
cat("size on disk:")
path <- file.path("data", "loans.rda")
dir.create("data", showWarnings = FALSE)
save(loans, file = path, compress = "xz")
cat("", format(file.size(path) / 1024, digits = 3), "KB\n")
