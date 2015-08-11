


class LongestCommonSubstringFinder:
"A class for finding the longest common substring of two or more strings"
    
    def findLongestCommonSubstring(self, string1, string2):
        """Returns the longest common substring of two parameter strings. If multiple common substrings
are of equal length and are the longest common substrings there as no guarantee as to which
will be returned"""
        countingArray = [[0]*len(string2) for i in range(len(string1))]
        ##we set up an array to locate the possible substrings
        lcss = ""
        for row in range(len(countingArray)):
            for column in range(len(countingArray[row])):
                if string1[row] == string2[column]:
                    countingArray[row][column] += 1 ##if there's a match the common string length is at least 1
                    if not (row == 0 or column == 0):
                        countingArray[row][column] += countingArray[row-1][column-1] ##memoize the length of the common string so far
                    if countingArray[row][column] > longestLength: ##if it's longer than our current longest string, we'll plan on returning it
                        longestLength = countingArray[row][column]
                        lcss = string1[row-longestLength+1:row+1]
        return lcss



    def findAllCommonSubstrings(self, string1, string2):
        """"Returns all common substring of two parameter strings."""
        countingArray = [[0]*len(string1) for i in range(len(string2))]
        ##we set up an array to locate the possible substrings
        candidates = set()
        for row in range(len(countingArray)):
            for column in range(len(countingArray[row])):
                if string2[row] == string1[column]:
                    countingArray[row][column] += 1
                    if not (row == 0 or column == 0):
                        countingArray[row][column] += countingArray[row-1][column-1]
                    candidates.add(string2[row-countingArray[row][column]+1:row+1])
        candidates = list(candidates)
        candidates.sort(key=len)
        candidates = list(reversed(candidates))
        return candidates



    def findLongestCommonSubstringManyStrings(listOfStrings):
        """Returns the longest common substring of two or more parameter strings. If multiple common substrings
are of equal length and are the longest common substrings there as no guarantee as to which
will be returned"""
##find all the substrings of our base case. Theoretically we could narrow the list of
        ##candidate strings down further by only calling findAllCommonSubstrings on the
        ##first two strings in the list of strings, however that wouldn't change the
        ##asymptotic growth by any significant amount so we're just taking a brute force approach
    candidateSubstrings = set()
    for i in range(len(listOfStrings[0])):
        j=0
        k=j+i
        while k <=len(listOfStrings[0]):
            candidateSubstrings.add(listOfStrings[0][j:k])
            j=j+1
            k=k+1

    #convert to a list in descending order
    candidateSubstrings = list(candidateSubstrings)
    candidateSubstrings = sorted(candidateSubstrings, key=len, reverse=True)
    
    #now check the substring against all the other strings
    for substring in candidateSubstrings:
        for sequence in listOfStrings:
            if substring not in sequence:
                found = False
                break #as soon as the substring isn't in one sequence we can move onto the next substring
            found = True 
        if found:
            return substring
    return "" #if there are no common substrings, implicitly the empty string is the longest substring

##It's worth noting that this solution is in quadratic time, but this problem can be solved in linear time
##using ukkonen's algorithm.
      
