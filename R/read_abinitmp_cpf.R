#' Read common/combined log file into a tibble
#'
#' This is a fairly standard format for log files - it uses both quotes
#' and square brackets for quoting, and there may be literal quotes embedded
#' in a quoted string. The dash, "-", is used for missing values.
#'
#' @inheritParams read_delim
#' @export
#' @import stringr
#' @import tibble
#' @examples
#' read_log(readr_example("example.log"))
read_abinitmp_cpf <- function(file, skip = 0, n_max = Inf, progress = show_progress()) {

  col_type <- do.call(cols, tribble(
    ~VERSION   ,   ~N_ATOM   ,   ~N_FRAG  ,  ~ATOM_INFO ,   ~N_ELEC  , ~FRAG_CHARGE,
  #|-----------|-------------|------------|-------------|------------|-------------|
     ~BINIDNG  ,  ~DISTANCE  ,   ~DIPOLE  ,   ~BASIS    ,   ~STATE   ,  ~METHOD    ,
  #|-------------------------|--------------------------|--------------------------|
        ~PARAM_AO_POP_APRX   , ~PARAM_POINT_CAHRGE_APRX ,   ~PARAM_DIMER_ES_APRX   ,
  #|-------------------------|--------------------------|--------------------------|
    ~NUC_ENERGY ,~ELEC_ENERGY ,~TOTAL_ENREGY ,~MONOMER  ,   ~IFIE    , ~MULTI_BODY ,
  #|-----------|-------------|------------|-------------|------------|-------------|
        "c"    ,    "n"      ,     "n"    ,    "c"      ,     "n"    ,    "c"      ,
        "c"    ,    "c"      ,     "c"    ,    "c"      ,     "c"    ,    "c"      ,
              "n"            ,           "n"            ,           "n"            ,
        "n"    ,    "n"      ,     "n"    ,    "c"      ,     "c"    ,    "c"      ,

  ))

  res <- read_delimited(file, tokenizer_abinitmp_cpf(),
                        col_names = names(col_type$cols),
                        col_types = col_type,
                        skip = skip,
                        n_max = n_max,
                        progress = progress)

  if(nrow(res) == 0) return(res)

  res$FRAG_CHARGE <- list(str_squish(res$FRAG_CHARGE) %>% str_split(" ") %>% nth(1) %>% as.integer())
  res$ATOM_INFO <- list(res$ATOM_INFO %>% read_table(col_names = c("id", "symbol", "atom_type", "residue_name", "residue_id", "fragment_id", "x", "y", "z", "c_1", "c_2", "c_3", "c_4","c_5","c_6","c_7")))
  res$BINIDNG <- list(res$BINIDNG %>% read_fwf(col_positions = fwf_widths(c(5, 5), c("from", "to"))))
  res$DISTANCE <- list(res$DISTANCE %>% read_delim(" ", col_names = c("from", "to", "distance"), col_types = cols(from = "i", to = "i", distance = "n") ,trim_ws = TRUE))
  res$DIPOLE <- list(res$DIPOLE %>% read_delim(" ", col_names = c("col_1", "col_2","col_3","col_4","col_5","col_6") ,trim_ws = TRUE))
  res$MONOMER <- list(res$MONOMER %>% read_delim(" ", col_names = FALSE ,trim_ws = TRUE))

  param_ifie <- do.call(fwf_positions, tribble(
 ~start,~end,   ~col_names       ,
#------|----|--------------------,
      1,  24, "nuc"              ,
     25,  48, "hf_elec"          ,
     49,  72, "hf_stat"          ,
     73,  96, "mp2_ifie"         ,
     97, 120, "scs_mp2_ifie"     ,
    121, 144, "mp3_ifie"         ,
    145, 168, "scs_mp3_ifie"     ,
    169, 192, "hf_ifie_bsse"     ,
    193, 216, "mp2_ifie_bsse"    ,
    217, 240, "scs_mp_ifie_bsse" ,
    241, 264, "mp3_ifie_bsse"    ,
    265, 288, "scs_mp3_ifie_bsse",))

  type_ifie <- cols(
    nuc = "n",
    hf_elec = "n",
    hf_stat = "n",
    mp2_ifie = "n",
    scs_mp2_ifie = "n",
    mp3_ifie = "n",
    scs_mp3_ifie = "n",
    hf_ifie_bsse = "n",
    mp2_ifie_bsse = "n",
    scs_mp_ifie_bsse = "n",
    mp3_ifie_bsse = "n",
    scs_mp3_ifie_bsse = "n"
  )

  res$IFIE <- list(res$IFIE %>% read_fwf(param_ifie, type_ifie))

  res
}



