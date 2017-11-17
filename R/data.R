#' Responses to questionaires from patients affected by breast cancer.
#'
#' A dataset containing the responses of 121 patients to 30 questions
#' about their quality of life.
#'
#' @format A dataframe with 121 lines and 32 columns. A line represents
#' a patient and a column are information about the patients
#' \describe{
#'   \item{Id}{patient Id}
#'   \item{q1-q28}{responses to 28 questions with number of levels equals to 4}
#'   \item{q29-q30}{responses to 22 questions with number of levels equals to 7}
#' }
#' @source
#' The table was determined based on data associated with the package
#' available on: \url{https://cran.r-project.org/web/packages/QoLR/index.html}
"dataqol"


#' Informaton from patients affected by breast cancer.
#'
#' A dataset containing the responses of 40 patients to 30 questions
#' about their quality of life. Furthermore, a variable indicates if the patient
#' survived from the disease.
#'
#' @format A dataframe with 40 lines and 32 columns. A line represents
#' a patient and a column are information about the patients
#' \describe{
#'   \item{Id}{patient Id}
#'   \item{q1-q28}{responses to 28 questions with number of levels equals to 4}
#'   \item{q29-q30}{responses to 22 questions with number of levels equals to 7}
#'	 \item{death}{if the patient survived (1) or not (2)}
#' }
#' @source
#' The table was determined based on data associated with the package
#' available on: \url{https://cran.r-project.org/web/packages/QoLR/index.html}
"dataqol.classif"