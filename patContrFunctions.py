class pat:
    import numpy as np
    start=0
    end=0
    Nbits=64
    Val=0 
    maskNbits=pow(2,Nbits)-1
    LineN=0
    File=''
    ioctrl=0
    clkctrl=0
    
class pat_default:
    start=0
    end=0
    Nbits=64
    Val=0 
    maskNbits=pow(2,Nbits)-1
    LineN=0
    File=''
    ioctrl=0
    clkctrl=0
    out1=0
    out2=0
    
####FORMATTING FUNCTIONS########################################
def hexFormat(val,fill):
    v=hex(val)  #hexadecimal value
    v=(v.lstrip('-0x')).rstrip('L')  #remove leading 0x and - if there
    v=v.zfill(fill) #inserts zeros at the beginning
    return v
def binFormat(val,fill):
    v=bin(val)  #binary value
    v=(v.lstrip('-0b')).rstrip('L')   #remove leading 0x and - if there and trailing L
    v=v.zfill(fill) #inserts zeros at the beginning
    return v
def decFormat(val,fill):
    v=str(val)  #decimal value
    v=v.zfill(fill) #inserts zeros at the beginning
    return v
################################################################
  
################################################################
##PATTERN CONTROL FUNCTIONS##
def setbit(bit,word):
    maskBit=1<<bit
    word=word|maskBit
    return(word)

def clearbit(bit,word):
    maskBit=1<<bit
    maskBitN=pat.maskNbits-maskBit
    word=word&maskBitN
    return(word)

def SB(bit):
    pat.Val=setbit(bit,pat.Val)
    
def CB(bit):
    pat.Val=clearbit(bit,pat.Val)
 
def CBs(*args):
    for i in args:
        pat.Val=clearbit(i,pat.Val)

def SBs(*args):
    for i in args:
        pat.Val=setbit(i,pat.Val)
         
def setoutput(bit):
    pat.ioctrl=setbit(bit,pat.ioctrl)

def setinput(bit):
    pat.ioctrl=clearbit(bit,pat.ioctrl)

def setclk(bit):
    pat.clkctrl=setbit(bit,pat.clkctrl)

def setinputs(*args):
    for i in args:
        setinput(i)

def setoutputs(*args):
    for i in args:
        setoutput(i)
        
def setclks(*args):
    for i in args:
        setclk(i)

def pw():
    address=pat.LineN
    value=pat.Val
    paw='patword'
    patline=paw+' '+hexFormat(address,4)+' '+hexFormat(value,16)+'\n'
#    print(decFormat(address,3)+'   '+binFormat(value,10))
    pat.File+=patline
    pat.LineN+=1
#    print(patline)

def saveToFile(fname):
    patEnd1='patloop0 0400 0400\npatnloop0 0\npatloop1 0400 0400\npatnloop1 0\npatloop2 0400 0400\n'
    patEnd2='patnloop2 0\npatwait0 0400\npatwaittime0 0\npatwait1 0400\npatwaittime1 0\npatwait2 0400\npatwaittime2 0'
    f=open(fname,'w')
    ioct='patioctrl '+hexFormat(pat.ioctrl,16)+'\n'
    clkct='patclkctrl '+hexFormat(pat.clkctrl,16)+'\n'
    plims='patlimits '+hexFormat(pat.start,4)+' '+hexFormat(pat.LineN-1,4)+'\n'
    pF=pat.File+ioct+clkct+plims+patEnd1+patEnd2
    f.write(pF)
    f.close()

def CLOCKS(bit,times):
    for i in range(0,times):
        SB(bit);pw()
        CB(bit);pw()

def serializer(value,serInBit,clkBit,nbits,msbfirst=1):
    """serializer(value,serInBit,clkBit,nbits,msbfirst=1)
    Produces the .pat file needed to serialize a word into a shift register.
    value: value to be serialized
    serInBit: control bit corresponding to serial in 
    clkBit: control bit corresponding to the clock 
    nbits: number of bits of the target register to load
    msbfirst: if 1 pushes in the MSB first (default), 
              if 0 pushes in the LSB first
    It produces no output because it modifies directly the members of the class pat via SB and CB"""
    c=value
    CBs(serInBit,clkBit)
    pw() #generate intial line with clk and serIn to 0
    start=0;stop=nbits;step=1
    if msbfirst:
        start=nbits-1;stop=-1;step=-1 #reverts loop if msb has to be pushed in first
    for i in range(start,stop,step):
        if c & (1<<i): 
            SB(serInBit)
            pw()
        else:
            CB(serInBit)
            pw()
        SB(clkBit)
        pw()
        CB(clkBit)
        pw() 
    CBs(serInBit,clkBit)
    pw() #generate final line with clk and serIn to 0     



##END PATTERN CONTROL FUNCTIONS##
##################################################################
##OTHER STUFF##

def plotDist(v,num_bins):
    import numpy as np
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    # the histogram of the data
    n, bins, patches = ax.hist(v, num_bins, normed=1)
    # add a 'best fit' line
    mu=np.average(v)
    print('Average:'+str(mu)+'\n')
    sigma=np.std(v)
    print('Standard deviation:'+str(sigma)+'\n')
    y = mlab.normpdf(bins, mu, sigma)
    ax.plot(bins, y, '--')
    ax.set_xlabel('Smarts')
    ax.set_ylabel('Probability density')
    ax.set_title(r'Histogram of single value')
    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    plt.show()
#   fig.savefig("test.png")


def pushbits(nbits,c):
  code=''
  for i in range(0,nbits-1):
    if c& (1<<i): 
      code = '1' + code
      #print("1")
      #print("code "+str(code))
    else:
      code= '0'+ code
      #print("0")
      #print("code "+str(code))
  return(code)
####################################
## DAC PERFORMANCE EVALUATION FUNCTIONS
def diffNlin(x,y,refN,refP,nbits):
    import numpy as np
    import matplotlib.pyplot as plt
    print("###########################################")
    print("#### DIFFERENTIAL NONLINEARITY ANALYSYS ###")
    print("###########################################")
    nPoints=2**nbits
    Vs=(refP-refN)/(nPoints-1)   
    DiffLin=((y[1:nPoints]-y[0:nPoints-1])-Vs)
    print("Maximum differential nonlinearity value:"+str(max(DiffLin))+' V'+ " (="+str(max(DiffLin)/Vs)+' LSB)' )
    print("Minimumdifferential nonlinearity value:"+str(min(DiffLin))+' V'+ " (="+str(min(DiffLin)/Vs)+' LSB)' )    
    plt.figure(1)
    plt.subplot(211)
    plt.plot(DiffLin)
    plt.ylabel("Diff. nonlinearity (V)")
#    plt.xlabel("DAC input code")
#    plt.figure(2)
    plt.subplot(212)
    plt.plot(DiffLin/Vs,'r')
    plt.ylabel("Diff. nonlinearity (LSB)")
    plt.xlabel("DAC input code")
    plt.show()
    plt.close(1)


def absNlin(x,y,refN,refP,nbits):
    import numpy as np
    import matplotlib.pyplot as plt
    nPoints=2**nbits
    Vs=(refP-refN)/(nPoints-1)
    gain=(refP-refN)/nPoints
    AbsLinRef=np.linspace(refN,refP,nPoints)
    DiffALR=AbsLinRef-y
    print("Maximum integral nonlinearity value:"+str(max(DiffALR))+' V'+ " (="+str(max(DiffALR)/Vs)+' LSB)' )
    print("Minimum integral nonlinearity value:"+str(min(DiffALR))+' V'+ " (="+str(min(DiffALR)/Vs)+' LSB)' )
    plt.figure(1)
    plt.subplot(211)
    plt.ylabel("Integ. nonlinearity(V)")
    plt.plot(DiffALR)
    plt.subplot(212)
    plt.ylabel("Integ. nonlinearity (LSB)")
#    plt.figure(2)
#    plt.xlabel("DAC input code")
#    plt.ylabel("Absolute linearity (LSB)")
    plt.xlabel("DAC input code")
    plt.plot(DiffALR/Vs,'r')
#    plt.plot(DiffALR/y)
    zeroError=y[0]-refN
    gainError=gain-(y[nPoints-1]-y[0])/nPoints
    fullScaleError=y[nPoints-1]-refP
    print("Zero error: "+str(zeroError))
    print("Full scale error: "+str(fullScaleError))
    print("Gain error: "+str(gainError))
    gainRelError=gainError/gain
    print("Gain relative error: "+str(gainRelError))
    plt.show()
    plt.close(1)

def DACperrformance(x,y,refN,refP,nbits):
    import numpy as np
    import matplotlib.pyplot as plt
    nPoints=2**nbits
    Vs=(refP-refN)/(nPoints-1)
    gain=(refP-refN)/nPoints
    AbsLinRef=np.linspace(refN,refP,nPoints)
    DiffALR=AbsLinRef-y
    print("Maximum integral Absolute nonlinearity value:"+str(max(DiffALR))+' V'+ " (="+str(max(DiffALR)/Vs)+' LSB)' )
    print("Minimum integral Absolute nonlinearity value:"+str(min(DiffALR))+' V'+ " (="+str(min(DiffALR)/Vs)+' LSB)' )
    plt.figure(1)
    plt.subplot(211)
    plt.ylabel("Integ. Absol. nonlinearity(V)")
    plt.plot(DiffALR)
    plt.subplot(212)
    plt.ylabel("Integ. Absol. nonlinearity (LSB)")
#    plt.figure(2)
#    plt.xlabel("DAC input code")
#    plt.ylabel("Absolute linearity (LSB)")
    plt.xlabel("DAC input code")
    plt.plot(DiffALR/Vs,'r')
#    plt.plot(DiffALR/y)
    zeroError=y[0]-refN
    gainError=gain-(y[nPoints-1]-y[0])/nPoints
    fullScaleError=y[nPoints-1]-refP
    print("Zero error: "+str(zeroError))
    print("Full scale error: "+str(fullScaleError))
    print("Gain error: "+str(gainError))
    gainRelError=gainError/gain
    print("Gain relative error: "+str(gainRelError))
    diffNlin(x,y,refN,refP,nbits)
    



