
def gen():
  k1=.1
  count=0
  f=open("parameters.txt",'w')
  while k1 < 11:
    m1=.001
    while m1<101:
      k2=.1
      while k2<11:
	m2=.001
	while m2 <101:
	  k3=.1
	  while k3 < 11:
	    m3=.001
	    while m3 < 101:
	      k4=.1
	      while k4 <11:
		m4=.001
		while m4 <101:
		  k5=.1
		  while k5 < 11:
		    m5=.001
		    while m5 <101:
		      k6=.1
		      while k6<11:
			m6=.001
			while m6<101:
			  x=str([k1,k2,k3,k4,k5,k6,m1,m2,m3,m4,m5,m6])
			  f.write(x)			  
			  f.write('\n')
			  m6=m6*10
			k6=k6*10
		      m5=m5*10		    
		    k5=k5*10
		  m4=m4*10
		k4=k4*10
	      m3=m3*10
	    k3=k3*10
	  m2=m2*10
	k2=k2*10
      m1=m1*10
    k1=k1*10
    
    print(count)
  f.close()

   
