//
//  NF_Child.java
//  LVTOOLS
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	Object NF_Child
//

//package NUCFAMTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class NF_Child implements Comparable{

	//Identity
	String Name;
	int Sex;			//1=male, 2=female
	NF_NucFam Family;	//this is a child in Family
	int nbSibs;
	int maxSibs;
	String[] sibnames;

	//Data
	String HapInfoMother;	//where father is homozygous
	String HapInfoFather;	//where mother is homozygous
	
	//Results
	int[][] RecodeF;
	int[][] RecodeM;
	boolean recode;
	
	int[] CrossoversF;
	int[] CrossoversM;
	boolean errors;

	TreeSet Errors;
	
	//Changes CANCER
	int[][] RecodeFerror ;
	int[][] RecodeMerror ;
	
	
	boolean chrx = false;
	
	
	//Constructor
	public NF_Child(String name, String data, int sex, int TotChildren, NF_NucFam famP, boolean x)
	{
	
		chrx = x;
		
		Name = name;
		Sex = sex;
		Family = famP;
		nbSibs = TotChildren - 1;
		sibnames = new String[nbSibs];
		maxSibs = nbSibs/2+1;
		
		HapInfoFather = "";
		HapInfoMother = "";
		
		SetHapInfo(data);
		
		RecodeF = null;
		RecodeM = null;
		CrossoversF = null;
		CrossoversM = null;
		
		Errors = null;
		
		//Changes CANCER
		RecodeFerror = null;
		RecodeMerror = null;
		
		recode = false;
		errors = false;
	}
	
	
	
	//Set infomrative haplotypes according to mother and father --> This child only is needed
	private void SetHapInfo(String data)
	{
		
		//For each informative marker, the child has two alleles. One is the mothers allele (homozygous mother) the other one is:
		// 1- the same allele if it is not the one on Father.HapInformative (child is homozygous)
		// 2- the fathers allele if it is the allele on Father.HapInformative (child is heterzygous)
		
		String FInfo = Family.Father.HapInformative;
		String MInfo = Family.Mother.HapInformative;
		
		String[] Data = data.split("\t");
		
		if(Data.length != FInfo.length()) Family.exit("Problem in SetHapInfo() with Data Format for child "+Name+" in Nuclear Family "+Family.ID);
		if(Data.length != MInfo.length()) Family.exit("Problem in SetHapInfo() with Data Format for child "+Name+" in Nuclear Family "+Family.ID);
		
		for(int i=0; i<Data.length;i++)
		{
			//Father
			if(FInfo.charAt(i)!='r')
			{
				String[] nuc = Data[i].split(" ");
				if(nuc[0].charAt(0)==nuc[1].charAt(0))// homozygote child
				{
					HapInfoFather += nuc[0].charAt(0);
				}
				else
				{
					HapInfoFather += FInfo.charAt(i);
				}
			}
			else//not informative
			{
				HapInfoFather += 'r';
			}
			
			//Mother
			if(MInfo.charAt(i)!='r')
			{
				String[] nuc = Data[i].split(" ");
				if(nuc[0].charAt(0)==nuc[1].charAt(0))
				{
					HapInfoMother += nuc[0].charAt(0);
				}
				else
				{
					HapInfoMother += MInfo.charAt(i);
				}
			}
			else
			{
				HapInfoMother += 'r';
			}
		}
		
	}
	
	
	//Subtable from Recode[][f] to Recode[][l] included
	private int[][] smallRecode(int[][] Recode, int f, int l)
	{
		int[][] newTable = new int[Recode.length][l-f+1];
		for(int i=0; i<newTable.length; i++)
		{
			for(int j=f; j<=l; j++)
			{
				newTable[i][j-f] = Recode[i][j];
			}
		}
		return newTable;
	}
	
	
	//Recode the siblings according to informative haplotype --> Siblings have to be set before recoding
	public void Recoding()
	{
		//Father
		
		RecodeF = new int[nbSibs][HapInfoFather.length()];
		
		int countSib = 0;
		if(!chrx)
		{
		
		//each child
		
		for (Iterator itr = Family.Children.iterator(); itr.hasNext();)
		{
			NF_Child act = (NF_Child)itr.next();
			if(!act.Name.equals(Name))
			{
				RecodeF[countSib] = RecodeIndiv( HapInfoFather, act.HapInfoFather);
				sibnames[countSib] = act.Name; // only once, no need to do this in mother
				countSib++;
			}

		}
		

		//when one missing (has 0 when != r), set everyone to missing.
		for(int i=0; i<RecodeF[0].length; i++)
		{
			if(HapInfoFather.charAt(i) != 'r') //informative
			{
				for(int j=0; j<RecodeF.length; j++)
				{
					if(RecodeF[j][i]==0)
					{
						for(int k=0; k<RecodeF.length; k++){
							RecodeF[k][i]=0;
						}
						j=RecodeF.length;
					}
				}
			}
		}
			
		}
		
		
		//Mother
		
		RecodeM = new int[nbSibs][HapInfoMother.length()];
		
		//each child
		countSib = 0;
		for (Iterator itr = Family.Children.iterator(); itr.hasNext();)
		{
			NF_Child act = (NF_Child)itr.next();
			if(!act.Name.equals(Name))
			{
				RecodeM[countSib] = RecodeIndiv( HapInfoMother, act.HapInfoMother);
				countSib++;
			}

		}
		
		//when one missing (has 0 when != r), set everyone to missing.	//DO I REALLY NEED THIS?? yes if one child missing but parent genotyped
		for(int i=0; i<RecodeM[0].length; i++)
		{
			if(HapInfoMother.charAt(i) != 'r')
			{
				for(int j=0; j<RecodeM.length; j++)
				{
					if(RecodeM[j][i]==0)
					{
						for(int k=0; k<RecodeM.length; k++){
							RecodeM[k][i]=0;
						}
						j=RecodeM.length;
					}
				}
			}
		}
		
		recode = true;
	
	}
	
	private int[] RecodeIndiv(String reference, String recodes)
	{
		if(reference.length() != recodes.length()) NF_NucFam.exit("Can't recode : reference sequence and recode sequence are not of same length");
		
		int[] recodeT = new int[recodes.length()];
		
		for(int i=0; i<recodes.length(); i++)
		{
			if(recodes.charAt(i)=='0' || reference.charAt(i)=='0' )
			{
				recodeT[i] = 0;
			}
			else if(reference.charAt(i)=='r' )
			{
				recodeT[i] = 0;
			}
			else if(reference.charAt(i)==recodes.charAt(i))
			{
				recodeT[i] = 1;
			}
			else
			{
				recodeT[i] = 2;
			}
		}
		
		return recodeT;
		
	}
	
	
	public void printRecode(char parent)
	{
		
		
		if (parent=='F')
		{
			if(RecodeF == null) NF_NucFam.exit("Undefined recoding : RecodeF == null");
			
			System.err.println("Recode table for child "+Name+" according to father data." );
			for(int i=0; i<RecodeF.length; i++)
			{
				for(int j=0; j<RecodeF[0].length; j++)
				{
					System.err.print(RecodeF[i][j]+" ");
				}
				System.err.println();
			}
		
		}
		
		else if (parent=='M')
		{
			if(RecodeM == null) NF_NucFam.exit("Undefined recoding : RecodeM == null");
			
			System.err.println("Recode table for child "+Name+" according to mother data." );
			for(int i=0; i<RecodeM.length; i++)
			{
				for(int j=0; j<RecodeM[0].length; j++)
				{
					System.err.print(RecodeM[i][j]+" ");
				}
				System.err.println();
			}
		
		}
		
		else{
			NF_NucFam.exit("Invalid entry "+parent +" in method printRecode()");
		}
	
	}
	
	public void printRecodeError(char parent)
	{
		if(Family.AllErrors == null) //meaning that FindErrors was not called previously
		{
			System.err.println("Error Set in family "+Family.ID+" is not up to date. Method printRecode() is called");
			printRecode(parent);
		}
		
		//New Recoding table without AllErrors found in all children.
		//Does not change RecodeF and RecodeM (Crossovers() still can be called)
		//int[][] RecodeFerror = null;
		//if(!chrx) {RecodeFerror = RemoveErrors(RecodeF,HapInfoFather);}
		//int[][] RecodeMerror = RemoveErrors(RecodeM,HapInfoMother);
		
		//Changes CANCER
		if(!chrx) {RecodeFerror = RemoveErrors(RecodeF,HapInfoFather);}
		RecodeMerror = RemoveErrors(RecodeM,HapInfoMother);
		
		
		if (parent=='F')
		{
			if(RecodeFerror == null) NF_NucFam.exit("Undefined recoding : RecodeFerror == null");
			
			System.err.println("Recode table for child "+Name+" according to father data." );
			for(int i=0; i<RecodeFerror.length; i++)
			{
				for(int j=0; j<RecodeFerror[0].length; j++)
				{
					System.err.print(RecodeFerror[i][j]+" ");
				}
				System.err.println();
			}
			
		}
		
		else if (parent=='M')
		{
			if(RecodeMerror == null) NF_NucFam.exit("Undefined recoding : RecodeMerror == null");
			
			System.err.println("Recode table for child "+Name+" according to mother data." );
			for(int i=0; i<RecodeMerror.length; i++)
			{
				for(int j=0; j<RecodeMerror[0].length; j++)
				{
					System.err.print(RecodeMerror[i][j]+" ");
				}
				System.err.println();
			}
			
		}
		
		else{
			NF_NucFam.exit("Invalid entry "+parent +" in method printRecodeError()");
		}
		
	}
	
	public ArrayList printRecodeErrorArray(char parent)
	{
		
				
		if(Family.AllErrors == null) //meaning that FindErrors was not called previously
		{
			System.err.println("Error Set in family "+Family.ID+" is not up to date. Method printRecode() is called");
			printRecode(parent);
		}
		
		//New Recoding table without AllErrors found in all children.
		//Does not change RecodeF and RecodeM (Crossovers() still can be called)
		//int[][] RecodeFerror = null;
		//if(!chrx) {RecodeFerror = RemoveErrors(RecodeF,HapInfoFather);}
		//int[][] RecodeMerror = RemoveErrors(RecodeM,HapInfoMother);
		
		//Changes CANCER
		if(!chrx) {RecodeFerror = RemoveErrors(RecodeF,HapInfoFather);}
		RecodeMerror = RemoveErrors(RecodeM,HapInfoMother);
		
		
		ArrayList SibsSharing = new ArrayList(RecodeFerror.length);//paternal for all siblings or maternal for all siblings
		
		if (parent=='F')
		{
			if(RecodeFerror == null) NF_NucFam.exit("Undefined recoding : RecodeFerror == null");
			
			//System.err.println("Recode table for child "+Name+" according to father data." );
			for(int i=0; i<RecodeFerror.length; i++)
			{
				String a = "F "+ Name + " " +sibnames[i]+" ";
				for(int j=0; j<RecodeFerror[0].length; j++)
				{
					a+=RecodeFerror[i][j]+" ";
				}
				SibsSharing.add(a);
				//System.err.println();
			}
			
		}
		
		else if (parent=='M')
		{
			if(RecodeMerror == null) NF_NucFam.exit("Undefined recoding : RecodeMerror == null");
			
			//System.err.println("Recode table for child "+Name+" according to mother data." );
			for(int i=0; i<RecodeMerror.length; i++)
			{
				String a = "M "+ Name + " " +sibnames[i]+" ";
				for(int j=0; j<RecodeMerror[0].length; j++)
				{
					a+=RecodeMerror[i][j]+" ";
				}
				SibsSharing.add(a);
				//System.err.println();
			}
			
		}
		
		else{
			NF_NucFam.exit("Invalid entry "+parent +" in method printRecodeError()");
		}
		
		return SibsSharing;
		
	}
	
	
	
	//Counting crossovers based on the recoding --> Siblings have to be set and recoding done before finding crossovers
	//Does not take genotyping errors into account : call CrossoversErrors() to account for genotyping errors
	public void Crossovers()
	{
		if(!chrx)
		{
		//Father
			if(RecodeF == null) NF_NucFam.exit("Undefined recoding (RecodeF == null). Cannot count crossovers.");
			CrossoversF = countCrossovers(RecodeF,HapInfoFather);
		}
		else 
		{
			CrossoversF = null;
		}
		
		//Mother
		if(RecodeM == null) NF_NucFam.exit("Undefined recoding (RecodeM == null). Cannot count crossovers.");
		CrossoversM = countCrossovers(RecodeM,HapInfoMother);
		
		//here, errors are not taken into account to compute crossovers tables
		errors = false;
		
	}
	
	private int[] countCrossovers(int[][] Recode, String HapInfo)
	{
	
		int[] cross = new int[HapInfo.length()];
		int[] last = new int[Recode.length];
		
		for(int i=0; i<last.length; i++)
		{
			last[i] = -1;
		}
		
		for(int i=0; i<Recode[0].length; i++)
		{
			
			//if(HapInfo.charAt(i)!='r' && Recode[0][i]!=0)//this site is informative
			if(HapInfo.charAt(i)!='r' && Recode[0][i]>0) //change 03-12-2013 because errors are -2
			{
				int switches = 0;//= new int[Recode.length];
				
				for(int j=0; j<Recode.length; j++)
				{
					int actstate = Recode[j][i];
					int laststate = last[j];
					
					if(laststate == -1)	//first time : setting last[j];
					{
						last[j]=actstate;
					}
					else if(laststate == 0)	//impossible
					{
						NF_NucFam.exit("laststate == 0 in Ccrossovers() (should never happen!!)");
					}
					else if(laststate!=actstate)
					{
						switches++;
						last[j]=actstate;	//new last[j]
					}
				}
				
				if(switches==nbSibs)// || switches >= maxSibs) **** VERY CONSERVATIVE!
				{
					cross[i] = 1;
				}
			}
		}
		
		return cross;
	}
	
	private int countCtable(int[] crossTable)
	{
		int sum = 0;
		for(int i=0; i<crossTable.length; i++){
			sum = sum + crossTable[i];
		}
		return sum;
	}
	
	public String printCrossovers2(char parent, int Kprint)
	{
		
		if(Kprint==0) return printCrossovers(parent);
		
		ArrayList Print = new ArrayList();
		
		if (parent=='F')
		{
			if(CrossoversF == null) NF_NucFam.exit("Undefined crossovers : CrossoversF == null");
			
			Print.add("Crossover events for child "+Name+" in family "+Family.ID+" for paternal meiosis ("+ Family.Father.Name +") with K = "+Kprint+" :" );
			
			String temp = null;
			int lastpos = -1;
			int track = 0;
			
			for(int i=0; i<HapInfoFather.length(); i++)
			{
				if(HapInfoFather.charAt(i) != 'r' && RecodeF[0][i]!=0) //informative
				{
					if(lastpos == -1)
					{
						lastpos = Family.AllPositions[i];
						//first informative position : no crossing over here
						if(CrossoversF[i] == 1) NF_NucFam.exit(lastpos+" is the first informative position CrossoversF[lastpos] can't be 1 (should never happen)");
						track = Kprint; //track goes from 0 to Kprint
					}
					else{
						if(CrossoversF[i] == 1) 
						{
							if(track >= Kprint)
							{
								if(temp!=null)
								{
									Print.add(temp);
									temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									//Print.add("\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i]);
									track = 1;
								}
								else
								{
									temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									track = 1;
								}
							}
							else
							{
								temp = null;
								track = Kprint;
							}
						}
						else
						{
							track++;
						}
						lastpos = Family.AllPositions[i];
					}
				}
			}
		}
		
		else if (parent=='M')
		{
			if(CrossoversM == null) NF_NucFam.exit("Undefined crossovers : CrossoversM == null");

			
			Print.add("Crossover events for child "+Name+" in family "+Family.ID+" for maternal meiosis ("+ Family.Mother.Name +")  with K = "+Kprint+" :" );
			
			String temp = null;
			int lastpos = -1;
			int track = 0;
			
			for(int i=0; i<HapInfoMother.length(); i++)
			{
				if(HapInfoMother.charAt(i) != 'r' && RecodeM[0][i]!=0) //informative
				{
					if(lastpos == -1)
					{
						lastpos = Family.AllPositions[i];
						//first informative position : no crossing over here
						if(CrossoversM[i] == 1) NF_NucFam.exit(lastpos+" is the first informative position CrossoversF[lastpos] can't be 1 (should never happen)");
						track = Kprint; //track goes from 0 to Kprint
					}
					else{
						if(CrossoversM[i] == 1) 
						{
							if(track >= Kprint)
							{
								if(temp!=null)
								{
									Print.add(temp);
									temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									//Print.add("\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i]);
									track = 1;
								}
								else
								{
									temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									track = 1;
								}
							}
							else
							{
								temp = null;
								track = Kprint;
							}
						}
						else
						{
							track++;
						}
						lastpos = Family.AllPositions[i];
					}
				}
			}
		
		}
		
		else{
			NF_NucFam.exit("Invalid entry "+parent +" in method printCrossovers()");
		}
		
		String myPrint = "";
		
		if(Print.size() > 1)
		{
			for(int i=0; i<Print.size(); i++)
			{
				if(i==0)
				{
					myPrint += (String)Print.get(i)+"\t"+(Print.size()-1) +"\n";
				}
				else
				{
					myPrint += (String)Print.get(i) + "\n";
				}
			}
			myPrint += "\n";
		}
		
		return myPrint;
	
	}
	
	
	
	public String printCrossovers(char parent, int Kprint)
	{
		
		//if(Kprint==0) return printCrossovers(parent);
		
		ArrayList Print = new ArrayList();
		
	
		ArrayList Blocs = new ArrayList(); //blocs of markers inbetween recombination events : this list will change depending on k
		
		
		if (parent=='F')
		{
			if(CrossoversF == null) NF_NucFam.exit("Undefined crossovers : CrossoversF == null");
			
			Print.add("Crossover events for child "+Name+" in family "+Family.ID+" for paternal meiosis ("+ Family.Father.Name +") with K = "+Kprint+" :" );
			
			int start = -1;
			int end = -1;
			int count = 0;
			
			for(int i=0; i<HapInfoFather.length(); i++)
			{
				
				if(RecodeF[0][i]!=0 && HapInfoFather.charAt(i) != 'r' && CrossoversF[i] != 1) //info but not rec
				{
					if(start == -1)
					{
						start = i;
					}
					count++;
					end = i;
				}
				else if(RecodeF[0][i]!=0 && HapInfoFather.charAt(i) != 'r' && CrossoversF[i] == 1) //info and rec
				{
					
					if(start == -1)
					{
						NF_NucFam.exit(Family.AllPositions[i]+" is the first position, CrossoversF can't be 1 (should never happen)");
					}
					int[] newt = new int[3];
					newt[0] = start;
					newt[1] = end;
					newt[2] = count;
					
					Blocs.add(newt);
					//reset
					start = i;
					count = 1;
				}
			}
			if(Blocs.size()!=0)		//addition 06/09/2011 (necessary bcs the last event was not printed!
			{
				//last event
				int[] newt = new int[3];
				newt[0] = start;
				newt[1] = end;
				newt[2] = count;
				Blocs.add(newt);
			}
			
		}
		else if (parent=='M')
		{
			if(CrossoversM == null) NF_NucFam.exit("Undefined crossovers : CrossoversM == null");
			
			Print.add("Crossover events for child "+Name+" in family "+Family.ID+" for maternal meiosis ("+ Family.Mother.Name +") with K = "+Kprint+" :" );
			
			int start = -1;
			int end = -1;
			int count = 0;
			
			for(int i=0; i<HapInfoMother.length(); i++)
			{
				
				if(RecodeM[0][i]!=0 && HapInfoMother.charAt(i) != 'r' && CrossoversM[i] != 1) //info but not rec
				{
					if(start == -1)
					{
						start = i;
					}
					count++;
					end = i;
				}
				else if(RecodeM[0][i]!=0 && HapInfoMother.charAt(i) != 'r' && CrossoversM[i] == 1) //info and rec
				{
					
					if(start == -1)
					{
						NF_NucFam.exit(Family.AllPositions[i]+" is the first position, CrossoversM can't be 1 (should never happen)");
					}
					int[] newt = new int[3];
					newt[0] = start;
					newt[1] = end;
					newt[2] = count;
					
					Blocs.add(newt);
					//reset
					start = i;
					count = 1;
				}
			}
			if(Blocs.size()!=0) //addition 06/09/2011 (necessary bcs the last event was not printed!
			{
				//last event
				int[] newt = new int[3];
				newt[0] = start;
				newt[1] = end;
				newt[2] = count;
				Blocs.add(newt);
			}
			
		}
		else{
			NF_NucFam.exit("Invalid entry "+parent +" in method printCrossovers()");
		}
				
		//System.err.println(parent+" "+((int[])Blocs.get(0))[0]+"-"+((int[])Blocs.get(0))[1]+"-"+((int[])Blocs.get(0))[2]);
		
		ArrayList BlocsR = new ArrayList();
		
		boolean print=false;
		//if(parent=='F' && Name.equals("383T") && Kprint==2) print=true;
		
		if(Kprint!=0)
		{
			BlocsR = ReduceBlocs(Blocs, Kprint, print);
			//second tour to make sure blocs below Kprint were not kept
			//if(Kprint>4) BlocsR = ReduceBlocs(BlocsR, Kprint, print);
		}
		else
		{
			BlocsR = Blocs;
		}
		
		//System.err.println(parent+" "+((int[])BlocsR.get(0))[0]+"-"+((int[])BlocsR.get(0))[1]+"-"+((int[])BlocsR.get(0))[2]);
		
		
		//Print
		Print.add("\tmarkBefore\tposition1\tposition2\tmarkAfter");
		for(int i =0; i<BlocsR.size()-1; i++)
		{
			int[] thisbloc = (int[])BlocsR.get(i);
			int[] nextbloc = (int[])BlocsR.get(i+1);
			
			Print.add("\t"+thisbloc[2]+"\t"+Family.AllPositions[thisbloc[1]]+"\t"+Family.AllPositions[nextbloc[0]]+"\t"+nextbloc[2]);
			
		}
		
		
		String myPrint = "";
		
		if(Print.size() > 1)
		{
			for(int i=0; i<Print.size(); i++)
			{
				if(i==0)
				{
					myPrint += (String)Print.get(i)+"\t"+(Print.size()-2) +"\n";
				}
				else
				{
					myPrint += (String)Print.get(i) + "\n";
				}
			}
			myPrint += "\n";
		}
		
		return myPrint;
		
	}
	
	private ArrayList ReduceBlocs(ArrayList Blocs, int Kprint, boolean print)
	{
		//if(print)System.err.println("start");
		ArrayList BlocsR = new ArrayList();
		
		//boolean continu = true;
		boolean skip = false;
		int[] rdouble = null;
		int c = 0;
		int i = 0;
		
		//if(print) System.err.print(Blocs.size());
		
		while(i<Blocs.size())
		{
			
			
			int[] t = (int[])Blocs.get(i);
			
			//if(print) System.err.print(i+" ");
			//t[0] = start
			//t[1] = end
			//t[2] = count
			
			if(t[2]>Kprint)
			{
				//if(print) System.err.println(" >");
				BlocsR.add(c,t);
				c++;
				i++;
				skip=false;
				rdouble=null;
			}
			else	// we are going to merge something
			{
				if(i==Blocs.size()-1)
				{
					//if(print) System.err.println("l ");
					int[] tlast = (int[])BlocsR.get(c-1);
					int[] tnew = new int[3];
					if(rdouble==null)
					{
						tnew[0]=tlast[0];
						tnew[1]=t[1];
						tnew[2]=t[2]+tlast[2];
						BlocsR.set(c-1,tnew);
						i++;
					}
					else //if last one was a rdouble, we have to merge it with t[2] as a bloc and c-1 is left untouched
					{
						tnew[0]=rdouble[0];
						tnew[1]=t[1];
						tnew[2]=t[2]+rdouble[2];
						BlocsR.add(c,tnew);
						i++;
					}
					rdouble=null;
					skip=false;
				
				}
				
				else 
				{
					
					int[] tnext = (int[])Blocs.get(i+1);
					if(tnext[2]<t[2] && !skip) // next one is smaller and we havent skip the last one
					{
						//we put this one in BlocsR, knowing it will be changed at the next step
						BlocsR.add(c,t);
						c++;
						i++;
						skip=true;
						rdouble=null; //if it was not null, it has been ignored
					}
					else if(tnext[2]==t[2] && !skip) // next one is same size and we havent skip the last one
					{
						//we skip both this one and the next
						i=i+2;
						
						if(rdouble==null)
						{
							rdouble = new int[3];
							rdouble[0] = t[0];
							rdouble[1] = tnext[1];
							rdouble[2] = t[2]*2;
						}
						//if the one before was also a double, increase the rdouble table
						else {
							rdouble[1] = tnext[1];
							rdouble[2] = rdouble[2]+t[2]*2;
						}
						
						
					}
					else if(tnext[2]>t[2] || skip) // next one is larger or we have skip the last one : must merge!
					{
						//if(print) System.err.println(" <");
						//if(print)System.err.println(i+" "+BlocsR.size());
						int[] tnew = new int[3];
						
						if(c==0)
						{
							tnew[0]=t[0];
							tnew[1]=tnext[1];
							tnew[2]=t[2]+tnext[2];
							BlocsR.add(c,tnew);
							i=i+2;
							c++;
						}
						else
						{
							int[] tlast = (int[])BlocsR.get(c-1);
							
							if(rdouble==null)
							{
								tnew[0]=tlast[0];
								tnew[1]=tnext[1];
								tnew[2]=t[2]+tlast[2]+tnext[2];
								BlocsR.set(c-1,tnew);
								i=i+2;
							}
							else //if last one was a rdouble, we have to merge it with t[2] as a bloc and c-1 is left untouched
							{
								tnew[0]=rdouble[0];
								tnew[1]=t[1];
								tnew[2]=t[2]+rdouble[2];
								BlocsR.add(c,tnew);
								c++;
								i++;
								rdouble=null;
							}
						}
						
						skip=false;
					}
				}
			
			//if(i>Blocs.size()-1) continu=false;
			}
		
		}
		
		
		
		//if(print)System.err.println("BlocsR.size()"+BlocsR.size());
		
		
		//if(print)System.err.println("end");
		return BlocsR;
	
	}
	
	
	
	
	public String printCrossovers(char parent)
	{
		ArrayList Print = new ArrayList();
		
		if (parent=='F')
		{
			if(CrossoversF == null) NF_NucFam.exit("Undefined crossovers : CrossoversF == null");
			
			Print.add("Crossover events for child "+Name+" in family "+Family.ID+" for paternal meiosis ("+ Family.Father.Name +") with K = 0 :" );
			
			//String temp = null;
			int lastpos = -1;
			//int track = 0;
			
			for(int i=0; i<HapInfoFather.length(); i++)
			{
				if(HapInfoFather.charAt(i) != 'r' && RecodeF[0][i]!=0) //informative
				{
					if(lastpos == -1)
					{
						lastpos = Family.AllPositions[i];
						//first informative position : no crossing over here
						if(CrossoversF[i] == 1) NF_NucFam.exit(lastpos+" is the first informative position CrossoversF[lastpos] can't be 1 (should never happen)");
						//track = Kprint; //track goes from 0 to Kprint
					}
					else{
						if(CrossoversF[i] == 1) 
						{
							//if(track >= Kprint)
							//{
								//if(temp!=null)
								//{
									//Print.add(temp);
									//temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									Print.add("\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i]);
									//track = 1;
								//}
								//else
								//{
									//temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									//track = 1;
								//}
							//}
							//else
							//{
							//	temp = null;
							//	track = Kprint;
							//}
						}
						else
						{
							//track++;
						}
						lastpos = Family.AllPositions[i];
					}
				}
			}
		}
		
		else if (parent=='M')
		{
			if(CrossoversM == null) NF_NucFam.exit("Undefined crossovers : CrossoversM == null");

			
			Print.add("Crossover events for child "+Name+" in family "+Family.ID+" for maternal meiosis ("+ Family.Mother.Name +")  with K = 0 :" );
			
			String temp = null;
			int lastpos = -1;
			int track = 0;
			
			for(int i=0; i<HapInfoMother.length(); i++)
			{
				if(HapInfoMother.charAt(i) != 'r' && RecodeM[0][i]!=0) //informative
				{
					if(lastpos == -1)
					{
						lastpos = Family.AllPositions[i];
						//first informative position : no crossing over here
						if(CrossoversM[i] == 1) NF_NucFam.exit(lastpos+" is the first informative position CrossoversF[lastpos] can't be 1 (should never happen)");
						//track = Kprint; //track goes from 0 to Kprint
					}
					else{
						if(CrossoversM[i] == 1) 
						{
							//if(track >= Kprint)
							//{
								//if(temp!=null)
								//{
									//Print.add(temp);
									//temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									Print.add("\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i]);
									//track = 1;
								//}
								//else
								//{
									//temp = "\tbetween informative positions\t"+lastpos+"\t"+Family.AllPositions[i];
									//track = 1;
								//}
							//}
							//else
							//{
							//	temp = null;
							//	track = Kprint;
							//}
						}
						//else
						//{
						//	track++;
						//}
						lastpos = Family.AllPositions[i];
					}
				}
			}
		
		}
		
		else{
			NF_NucFam.exit("Invalid entry "+parent +" in method printCrossovers()");
		}
		
		String myPrint = "";
		
		if(Print.size() > 1)
		{
			for(int i=0; i<Print.size(); i++)
			{
				if(i==0)
				{
					myPrint += (String)Print.get(i)+"\t"+(Print.size()-1) +"\n";
				}
				else
				{
					myPrint += (String)Print.get(i) + "\n";
				}
			}
			myPrint += "\n";
		}
		
		return myPrint;
	
	}
	
	
	//Counting crossovers with errors correction made before --> Siblings have to be set and recoding done. 
	//Uses RemoveErrors() method to find errors and calls the countCrossovers() method at the end.
	public void CrossoversErrors()
	{
		
		if(Family.AllErrors == null) //meaning that FindErrors was not called previously
		{
			System.err.println("Error Set in family "+Family.ID+" is not up to date. Method Crossovers() is called");
			Crossovers();
		}
		
		//New Recoding table without AllErrors found in all children.
		//Does not change RecodeF and RecodeM (Crossovers() still can be called)
		//int[][] RecodeFerror = null;
		//if(!chrx) {RecodeFerror = RemoveErrors(RecodeF,HapInfoFather);}
		//int[][] RecodeMerror = RemoveErrors(RecodeM,HapInfoMother);
		
		//Changes CANCER
		if(!chrx) {RecodeFerror = RemoveErrors(RecodeF,HapInfoFather);}
		RecodeMerror = RemoveErrors(RecodeM,HapInfoMother);
		
		
		//Crossovers are changed. If boolean errors = false, crossovers are computed without taking errors into accout
		//here errors = true, crossovers are computed and errors are taken into account.
		if(!chrx) 
		{
			CrossoversF = countCrossovers(RecodeFerror,HapInfoFather);
		}
		else 
		{
			CrossoversF = null;
		}
		CrossoversM = countCrossovers(RecodeMerror,HapInfoMother);
		
		errors = true;
		
		//if Crossovers() is called, CrossoversF and CrossoversM will be computed without errors and errors will be set to false
	}
	
	public String FindErrors(int interval)
	{
		Errors = new TreeSet();
		
		if(RecodeF == null) NF_NucFam.exit("Undefined recoding (RecodeF == null). Cannot find errors.");
		if(RecodeM == null) NF_NucFam.exit("Undefined recoding (RecodeM == null). Cannot find errors.");
		
		//Fill TreeSet Errors
		return FindErrors(RecodeF,HapInfoFather,interval) + FindErrors(RecodeM,HapInfoMother,interval);
			
	}
	
	private String FindErrors(int[][] Recode, String HapInfo, int interval)
	{
		
		String str = "";
		
		int f = -1;
		int a = -1;
		int l = -1;
		
		//initiation
		boolean init = false;
		int m = 0;
		while(!init)	//*** verify that we have at least 3 informatives markers
		{
			//System.err.println("f="+f+"a="+a+"l="+l+init);
			if(HapInfo.charAt(m)!='r' && Recode[0][m]>0) //chnages 03-12-2013 Recode[0][m]!=0
			{
				if(f == -1) f=m; 
				else if(a == -1) a=m;
				else if(l == -1) 
				{
					l=m;
					init = true;
				}
			}
			m++;
			if(m==HapInfo.length()) init = true;
		}
		
		
		
		//step
		boolean end = false;
		
		if(l==-1) end = true;
		
		while(!end)
		{
			int[][] smallRecode = smallRecode(Recode,f,l);
			String smallHapInfo = HapInfo.substring(f,l+1);
			
			
			int[] with = countCrossovers(smallRecode, smallHapInfo);
			//for(int i=0;i<with.length;i++) {System.err.print(with[i]+" ");}
			//System.err.println();
			
			int new_a = a - f;
			for(int i=0; i<smallRecode.length; i++)
			{
				smallRecode[i][new_a] = 0;
			}
			int[] without = countCrossovers(smallRecode, smallHapInfo);
			//for(int i=0;i<with.length;i++) {System.err.print(without[i]+" ");}
			//System.err.println();
			
			
			
			//*** Only THIS case is considered as an error for the moment
			if(countCtable(with) == 2 && countCtable(without) == 0)  // case 10201 or 20102
			{
				//verify distances
				int fpos = Family.AllPositions[f];
				int apos = Family.AllPositions[a];
				int lpos = Family.AllPositions[l];
				
				if(apos-fpos < interval && lpos-apos < interval) //mean of a hotspot of recombination every 50kb
				{
					if(Recomb.ValidErrors==null || !((Recomb.ValidErrors).contains(""+apos)))
					{
						Errors.add(new Integer(a));
						str += "Error detected in child "+Name+" at SNP # "+a+" ("+apos+", interval = "+(lpos-fpos)+").\n";
					}
					else
					{
						str += "Error ignored in child "+Name+" at SNP # "+a+" ( "+apos+" ).\n";
					}
				}
			}
			//else if(countCtable(with) == 0 && countCtable(without) == 0){} // case 10101 or 20202
			//else if(countCtable(with) == 1 && countCtable(without) == 1){} // case 10102 or 20101 or 20201 or 10202
			//else{
			//	NF_NucFam.exit("Impossible pattern of recombination for "+Name+" (should never happen) : "+countCtable(with)+" "+countCtable(without)+" "+f+"-"+a+"-"+l);
			//}
			
			//next
			f = a;
			a = l;
			
			boolean find = false;
			int k = l+1;
			while(!find && k<HapInfo.length())
			{
				if(HapInfo.charAt(k)!='r' && Recode[0][k]>0) //changes 03-12-2013 Recode[0][k]!=0
				{
					l=k;
					find = true;
				}
				k++;
			}
			if(!find) end = true;
			
		}
		
		return str;
	}
	
	private int[][] RemoveErrors(int[][] Recode, String HapInfo)
	{
		if(Family.AllErrors == null) NF_NucFam.exit("Error Set in family "+Family.ID+" is not up to date.");
		
		for (Iterator itr = Family.AllErrors.iterator(); itr.hasNext();)
		{
			int i = ((Integer)itr.next()).intValue();
			if(HapInfo.charAt(i) != 'r' && Recode[0][i]>0) //changes 03-12-2013 Recode[0][k]!=0
			{
				for(int j=0; j<Recode.length; j++)
				{
					//Recode[j][i] = 0;//old one
					Recode[j][i] = -2; //changes 03-12-2013
					
				}
			}
		}
		
		return Recode;
		
	}
	
	
	
	//Comparable
	public int compareTo(Object that)
	{
		
		if (this == that) return 0;
		
		if(that == null){
			throw new NullPointerException();
		}
		
		if(!(that instanceof NF_Child)){
			throw new ClassCastException();
		}
		else
		{
			return ((NF_Child)this).Name.compareTo(((NF_Child)that).Name);
		}
		
	}
	
}

