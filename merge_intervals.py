#!/usr/bin/python

import sys


def compare(s1,s2):
	if (s1[2]>s2[2]):
		return 1
	elif (s1[2]<s2[2]):
		return -1
	if (s1[0]>s2[0]):
		return 1
	elif (s1[0]<s2[0]):
		return -1
	return 0 

def usage():
	print sys.argv[0] + " forgiveness"

def collapse(intervals):
	if (len(intervals)==0):
		return intervals
	temp=[]
	start=intervals[0][0]
	end=intervals[0][1]
	call=intervals[0][3]
	chr=intervals[0][2]
	length=end-start+1
	for i in intervals[1:]:
		next_start=i[0]
		next_end=i[1]
		next_call=i[3]
		next_length=next_end-next_start+1
		next_chr=i[2]
		if chr==next_chr and next_call*call>0 and next_start<=end:
			length=length+(max(next_end-end+1,0))
			call=(call*length+next_call*next_length)/length
			if ((next_end-end+1)>0):
				end=next_end
		else:
			temp.append([start,end,chr,call])
			start=next_start
			end=next_end
			call=next_call
			length=next_end-next_start+1
			chr=next_chr
	temp.append([start,end,chr,call])
	return temp


def merge_close(forgiveness,intervals):
        if (len(intervals)==0):
                return intervals
        temp=[]
        start=intervals[0][0]
        end=intervals[0][1]
        call=intervals[0][3]
        chr=intervals[0][2]
        length=end-start+1
        for i in intervals[1:]:
                next_start=i[0]
                next_end=i[1]
                next_call=i[3]
                next_length=next_end-next_start+1
                next_chr=i[2]
                #assert(chr==next_chr)
		#assert(next_start>end)
		gap_length=next_start-end - 1
		allowance=(float(forgiveness)*(next_end-start))/100
		if (gap_length<=allowance) and chr==next_chr and next_call*call>0:
			#need to merge
			call=(call*length+next_call*next_length+gap_length*(call/2+next_call/2))/(length+gap_length+next_length)
			length=length+gap_length+next_length
			end=next_end
		else:
                        temp.append([start,end,chr,call])
                        start=next_start
                        end=next_end
                        call=next_call
                        length=next_end-next_start+1
			chr=next_chr
        temp.append([start,end,chr,call])
        return temp



if __name__=='__main__':
	if len(sys.argv)!=2:
		usage()
		sys.exit(1)
	
	forgiveness = int(sys.argv[1])

	calls=[]
	#gains=[]
	#loses=[]
	
	for line in sys.stdin.readlines():	
		line=line.split()
		if (len(line)==4):
			chr=line[0]
			interval_start=int(line[1])
			interval_end=int(line[2])
			call=float(line[3])
			calls.append([interval_start,interval_end,chr,call])
			#if call>0:
			#	gains.append([interval_start,interval_end,chr,call])
			#else:
			#	loses.append([interval_start,interval_end,chr,call])
		else:
			print "#WARNING! Skipping line " + "\""+ " ".join(line) + "\""

	#gains.sort()
	#gains=collapse(gains)
	#gains=merge_close(forgiveness,gains)
	#loses.sort()
	#loses=collapse(loses)
	#loses=merge_close(forgiveness,loses)
	
	#results=gains+loses
	#results.sort()

	calls.sort(cmp=compare)
	calls=collapse(calls)
	calls=merge_close(forgiveness,calls)

	for r in calls:
		print r[2],r[0],r[1],r[3]

	#take out overlaps
		

	#do  gains
	#diff=True
	#while diff:
	#	diff=False
		


