import cPickle as pickle
from os import getenv

l = pickle.load(open(getenv('SUBMIT_WORKDIR')+'/submission.pkl'))
s = l[-1]
print s.cluster_id
statii = s.query_status()
print 'Job summary:'
for k,v in statii.iteritems():
      print '\t %10s : %5i'%(k,len(v))
      for j in v:
            if 'MET' in j.name:
                  print j.name

