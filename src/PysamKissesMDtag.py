
import pysam
import string
import random

def parseMD (md):
    '''
    Input:
        md: MD string reported by read.opt("MD")
    Output:
        MD array in pysam syntax:
            D : Deletion : 2
            = : Equal    : 7
            X : Diff     : 8
    Modified from Michiel de Hoon (http://genome.gsc.riken.jp/osc/english/dataresource/)
    '''
    # md = md.replace("MD:Z:", "")
    result = []
    number = ""
    deletion = ""
    for c in md:
        if c in "0123456789":
            if deletion:
                ## deletion results
                deletion = (2, (len(deletion)) - 1)
                result.append(deletion)
                deletion = ""
                number = ""
            number += c
        elif c=='^':
            ## matched results previous
            number = (7, int(number))
            result.append(number)
            deletion = "^"
            number = ""
        elif deletion:
            deletion += c
        else:
            if number!="":
                ## matched results
                number = (7, int(number))
                result.append(number)
            ## Mismatch results
            c = (8, len(c))
            result.append(c)
            number = ""
    assert deletion==""
    if number!="":
        ## matched results
        number = (7, int(number))
        result.append(number)
    return result

def truncateMD (queryLen, md):
    '''
    Step-wise truncate MD tag, report MD given query length
    Input:
        queryLen: [int], length of query
        md: [list], parseMD(md) object
    Output:
        list[0]: [list] MD in length of query
        list[1]: [list] truncated MD list
    '''
    accumLen = 0
    queryMD = list()
    tmpMD = md
    for i in xrange(len(md)) :
        # print "block" + str(i)
        blockLen = md[i][1]
        blockTyp = md[i][0]
        accumLen = accumLen + blockLen
        # print accumLen
        if (queryLen < accumLen):
            md[i] = (blockTyp, (accumLen - queryLen))
            tmpMD = md[i: ]
            queryMD.append((blockTyp, (queryLen - (accumLen - blockLen))))
            # print queryMD
            # print tmpMD
            return([queryMD, tmpMD])

        elif (accumLen == queryLen) :
            queryMD.append(md[i])
            tmpMD = md[(i + 1): ]
            # print queryMD
            # print tmpMD
            return([queryMD, tmpMD])
        elif (queryLen > accumLen):
            queryMD.append(md[i])
            tmpMD = md[(i + 1) : ]
            # print queryMD
            # print tmpMD
        else :
            pass
    if (queryLen > accumLen) :
        return("Warnings!!! Query Length larger than MD string length")

def moreCigar (cigar, md):
    '''
    Goal: Absorbing MD tag information to modify cigar information
    Note: pysam's cigar syntax (list)
    Demo: (in string syntax)
    5M + 3T4 -> 3=1X2=
    Input:
        cigar: read.cigar() object
        md   : read.opt("MD") object
    '''

    md = parseMD(md)

    output = list() ## modified cigar list
    tmpMd = md     ## truncated md list

    for i in xrange(len(cigar)) :
        blockTyp = cigar[i][0]
        blockLen = cigar[i][1]

        if blockTyp in [0, 2]:
            # print "Match/Mismatch, & Deletion:"
            tmpTrunc = truncateMD(blockLen, tmpMd)
            output = output + tmpTrunc[0]
            tmpMd = tmpTrunc[1]

        else:
            output.append(cigar[i])
    return(output)

# Demo of alignments
# bamFile = "/Users/yunyan/Projects/TestData/bin.md.e.bam"
# aln = pysam.Samfile(bamFile, "rb")
#
# I = 0
# nums = random.sample(xrange(10000), 20 )
# for read in aln:
#     I += 1
#     if I in nums:
#         print
#         print "Read " + str(I) + ":"
#         print read.cigarstring
#         # print read
#         # print read.alen
#         md = read.opt("MD")
#         print md
#         # print "cigar array: "
#         cigar = read.cigar
#         print cigar
#         print moreCigar(cigar, md)

# testMD = "1A1^AT40"
# parsedMd = parseMD(testMD)
# print "parsed Md: "
# print(parsedMd)
#
# queryLen = 10
# print(str(queryLen) + " query Length")
# truncatedMd = truncateMD(queryLen, parsedMd)
# print "Output: md list of query length"
# print (truncatedMd[0])
# print "truncated md list"
# print (truncatedMd[1])

print "################"
print "Demo: "
demoMD = "31^AT3T14"
print "demo MD tag: "
print demoMD
print "demo CIGAR List: "
cigarList = [(0, 31), (2,2), (3, 94), (0, 18)]
print cigarList
print "######"
print "cigar list with match/mismatch information: "
moreCigarList = moreCigar(cigarList, demoMD)
print moreCigarList

# #
# print "#######"
# debugMD = "36T14"
# debugMD = parseMD(debugMD)
# queryLen = 33
# print "1. " + str(queryLen) + "region: "
# trucateMD1 = truncateMD(queryLen, debugMD)
# print trucateMD1[0]
# print trucateMD1[1]
#
# queryLen = 18
# debugMD = trucateMD1[1]
# print "2. " + str(queryLen) + "region: "
# trucateMD2 = truncateMD(queryLen, debugMD)
# print trucateMD2[0]
# print trucateMD2[1]
