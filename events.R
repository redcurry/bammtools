# Returns an 'events' data frame by reading from the specified path.
# If n is specified, it uses only the last n generations of the data.
ReadEvents <- function(path, n = 0)
{
  if (n > 0)
    ReadEventsNGenerations(path, n)
  else
    ReadEventsAllGenerations(path)
}


ReadEventsNGenerations <- function(path, n)
{
  LastGenerations(ReadEventsAllGenerations(path), n)
}


ReadEventsAllGenerations <- function(path)
{
  # na.strings='na' prevents interpreting 'NA' as NA,
  # which in the data means an event occurs on a tip node.
  read.csv(path, col.names = EventColumnNames(),
    stringsAsFactors = F, na.strings = 'na')
}


EventColumnNames <- function()
{
  c('Generation', 'Taxon1', 'Taxon2', 'Time',
    'LambdaInit', 'LambdaRate', 'MuInit', 'MuRate')
}


LastGenerations <- function(events, n)
{
  generation.first <- FirstOfLast(unique(events$Generation), n)
  events[events$Generation >= generation.first, ]
}


FirstOfLast <- function(x, n)
{
  head(tail(x, n), 1)
}
