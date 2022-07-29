HBD.overlap <- function(i, pop_hbd, start, end) {
  hbd_start <- pop_hbd[i,'start_dist']
  hbd_end <- pop_hbd[i, 'end_dist']
  if (hbd_start > start & hbd_start < end | hbd_end > start & hbd_start < end) {
    return(paste0(pop_hbd[i, 'famid'],':', pop_hbd[i, 'id']))
  }
}