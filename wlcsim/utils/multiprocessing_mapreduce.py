import collections
import itertools
import multiprocessing

class SimpleMapReduce(object):

    def __init__(self, map_func, reduce_func, num_workers=None,
                 pool_type='process'):
        """
        map_func

          Function to map inputs to intermediate data. Takes as argument
          one input source and returns a list of objects that can be
          unpacked as a (key, value) pair to be reduced.

        reduce_func

          Function to reduce partitioned version of intermediate data to
          final output. Takes a single argument item=(key,valuelist)
          where key is as produced by map_func and valuelist is the
          sequence of the values associated with that key, one value for
          each worker that produced that key.

        num_workers

          The number of workers to create in the pool. Defaults to the
          number of CPUs available on the current host.
        """
        self.map_func = map_func
        self.reduce_func = reduce_func
        if pool_type == 'thread':
            self.pool = multiprocessing.pool.ThreadPool(num_workers)
        elif pool_type == 'process':
            self.pool = multiprocessing.Pool(num_workers)
        else:
            raise NameError('No such multiprocessing pool type: ' +
                            pool_type)

    def partition(self, mapped_values):
        """Organize the mapped values by their key.
        Returns an unsorted sequence of tuples with a key and a sequence of values.
        """
        partitioned_data = collections.defaultdict(list)
        for key, value in mapped_values:
            partitioned_data[key].append(value)
        return partitioned_data.items()

    def __call__(self, inputs, chunksize=1):
        """Process the inputs through the map and reduce functions given.

        inputs
          An iterable containing the input data to be processed.

        chunksize=1
          The portion of the input data to hand to each worker.  This
          can be used to tune performance during the mapping phase.
        """
        # get list of responses from each worker, which should each be lists of
        # tuples they mapped
        map_responses = self.pool.map(self.map_func, inputs, chunksize=chunksize)
        # combine the list of tuples from all workers into a master list, then
        # partition it based on key
        partitioned_data = self.partition(itertools.chain(*map_responses))
        # pass the .items() iterator with the promised key,valuelist pairs
        # to the reduce_func
        reduced_values = self.pool.map(self.reduce_func, partitioned_data)
        return reduced_values
