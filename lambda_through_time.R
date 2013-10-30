library(ape)

# Source files to include
source('events.r')


MeanLambdaThroughTime <- function(tree.path, event.path, n, n.breaks)
{
  tree <- read.tree(tree.path)
  events <- ReadEvents(event.path, n)
  events$Node <- EventNode(tree, events)
  times <- seq(0, MaxTime(tree), length.out = n.breaks)

  mean.lambda.at.times <- data.frame()

  generations <- unique(events$Generation)
  for (generation in generations)
  {
    generation.events <- events[events$Generation == generation, ]
    generation.events$EventId <- seq(nrow(generation.events))
    edge.events <- EdgeEvents(tree, generation.events)
    mean.lambda.at.times <- rbind(mean.lambda.at.times,
      MeanEventLambdaAtTimes(times, edge.events, events))
  }

  names(mean.lambda.at.times) <- times
  mean.lambda.at.times
}


MaxTime <- function(tree)
{
  max(branching.times(tree))
}


EventNode <- function(tree, events)
{
  MrcaNodes(tree, events$Taxon1, events$Taxon2)
}


MrcaNodes <- function(tree, taxon1, taxon2)
{
  taxon.pairs <- cbind(taxon1, taxon2)
  apply(taxon.pairs, 1, MrcaNode, tree = tree)
}


MrcaNode <- function(tree, tip.names)
{
  if (tip.names[2] == 'NA')
    TipNodes(tree, tip.names[1])
  else
    getMRCA(tree, TipNodes(tree, tip.names))
}


TipNodes <- function(tree, tip.names)
{
  Index(tree$tip.label, tip.names)
}


EdgeEvents <- function(tree, events)
{
  root.node <- length(tree$tip.label) + 1
  root.descendants <- DirectDescendants(tree, root.node)

  edge1 <- c(root.node, root.descendants[1])
  edge2 <- c(root.node, root.descendants[2])

  node.times <- NodeTimes(tree)

  # Assumption: there should only be one event at time 0.0
  root.event.id <- events[events$Time == 0, ]$EventId

  edge.events <- data.frame()

  edge.events <- RecurseCreateEdgeEvents(tree$edge, edge1, node.times,
    events, root.event.id, edge.events)

  edge.events <- RecurseCreateEdgeEvents(tree$edge, edge2, node.times,
    events, root.event.id, edge.events)

  edge.events
}


DirectDescendants <- function(tree, node)
{
  if (IsTip(tree, node))
    return(c())
  tree$edge[Index(tree$edge[, 1], node), 2]
}


IsTip <- function(tree, node)
{
  node <= length(tree$tip.label)
}


NodeTimes <- function(tree)
{
  ReverseTimes(c(TipTimes(tree$tip.label), branching.times(tree)))
}


TipTimes <- function(tip.label)
{
  n <- length(tip.label)
  times <- rep(0.0, n)
  names(times) <- seq(n)
  times
}


ReverseTimes <- function(times)
{
  max(times) - times
}


RecurseCreateEdgeEvents <- function(tree.edges, edge, node.times,
  events, event.id, edge.events)
{
  edge.times <- EdgeTimes(edge, node.times)

  edge.events.on.edge <- EdgeEventsOnEdge(edge, edge.times, events, event.id)
  edge.events <- rbind(edge.events, edge.events.on.edge)

  last.event.id <- Last(edge.events.on.edge)$EventId

  edge.descendants <- EdgeDescendants(edge, tree.edges)
  if (nrow(edge.descendants) > 0) {
    for (row.id in seq(nrow(edge.descendants))) {
      edge.descendant <- edge.descendants[row.id, ]
      edge.events <- RecurseCreateEdgeEvents(tree.edges, edge.descendant,
        node.times, events, last.event.id, edge.events)
    }
  }

  edge.events
}


EdgeTimes <- function(edge, node.times)
{
  c(node.times[edge[1]], node.times[edge[2]])
}


EdgeEventsOnEdge <- function(edge, edge.times, events, starting.event.id)
{
  # TODO: update starting.event.id if there is an event right on start time

  events.on.edge <- EventsOnEdge(edge, edge.times, events)
  events.on.edge <- events.on.edge[order(events.on.edge$Time), ]

  event.time.begin <- edge.times[1]
  event.id <- starting.event.id

  edge.events <- data.frame()
  if (nrow(events.on.edge) > 0) {
    for (row.id in seq(nrow(events.on.edge))) {
      event.on.edge <- events.on.edge[row.id, ]
      edge.events <- rbind(edge.events,
        c(event.time.begin, event.on.edge$Time, event.id))
      event.time.begin <- event.on.edge$Time
      event.id <- event.on.edge$EventId
    }
  }

  edge.events <- rbind(edge.events,
    c(event.time.begin, edge.times[2], event.id))

  names(edge.events) <- c('EventTimeBegin', 'EventTimeEnd', 'EventId')
  edge.events
}


EventsOnEdge <- function(edge, edge.times, events)
{
  events[(events$Node == edge[2]) & (events$Time < edge.times[2]), ]
}


EdgeDescendants <- function(edge, tree.edges)
{
  tree.edges[tree.edges[, 1] == edge[2], ]
}


MeanEventLambdaAtTimes <- function(times, edge.events, events)
{
  sapply(times, MeanEventLambdaAtTime, edge.events, events)
}


MeanEventLambdaAtTime <- function(time, edge.events, events)
{
  edge.event.index <- EdgeEventAtTime(time, edge.events)
  mean(EventLambdaAtTime(time, edge.event.index, events))
}


EdgeEventAtTime <- function(time, edge.events)
{
  edge.events[EdgeEventCrossesTime(time, edge.events), ]$EventId
}


EdgeEventCrossesTime <- function(time, edge.events)
{
  (edge.events$EventTimeBegin <= time) & (time < edge.events$EventTimeEnd)
}


EventLambdaAtTime <- function(time, event.index, events)
{
  LambdaAtTime(time - events[event.index, ]$Time,
    events[event.index, ]$LambdaInit, events[event.index, ]$LambdaRate)
}


LambdaAtTime <- function(time, lambda.init, lambda.rate)
{
  lambda.init * exp(lambda.rate * time)
}


# Helper functions


Last <- function(x)
{
  tail(x, 1)
}


Index <- function(x, values)
{
  sapply(values, WhichIndex, x = x)
}


WhichIndex <- function(x, value)
{
  which(x == value)
}
