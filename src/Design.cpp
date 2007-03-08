#include "Design.h"

//public member function:

void Design::setArray(vector<int> a)
{
	Array = a; 
}

void Design::setCond(vector<int> c) 
{
	Cond = c;
}

void Design::setDye(vector <int> d)
{
	Dye = d;
}

vector <int> Design::getArray() 
{
	return Array;
}

vector <int> Design::getCond() 
{
	return Cond;
}

vector <int> Design::getDye() 
{
	return Dye;
}

vector <Design> Design::splitBlock(Design des) 
{
	vector <int> a = des.getArray();
	vector <int> c = des.getCond();
	vector <int> d = des.getDye(); 

	// find out unique conditions.

	vector <int> c_unique(c.begin(),c.end());
	sort(c_unique.begin(),c_unique.end());
	vector <int>::iterator pc_unique = unique(c_unique.begin(),c_unique.end());
	c_unique.erase(pc_unique,c_unique.end());

	vector <int>::iterator pa = a.begin();
	vector <int>::iterator pc = c.begin();
	vector <int>::iterator pd = d.begin();
	pc_unique = c_unique.begin();

	vector <Design> designblock;

	if (c_unique.size() > 1)
	{
		int currentcond = 0;
		int currentarray = 0;
		int currentdye = 0;
		
		unsigned int i = 0;
		unsigned int j = 0;

		while (pc_unique!=c_unique.end())
		{
			currentcond = *pc_unique;

			vector <int> tmpc;
			tmpc.push_back(currentcond);
					
			for (i=0;i<tmpc.size();i++)
			{
				currentcond = tmpc.at(i);

				int index = 0;
				vector <int> vindex;  // vector vindex is for storing corresponding index in which the element is currentcond.
				pc = c.begin();

				while (pc != c.end())  // if element in index vector == currentcond, store the index in vindex.
				{
					if (*pc == currentcond)
					{
						vindex.push_back(index);
					}

					index ++;
					pc ++;
				}

				vector <int>::iterator pindex = vindex.begin();

				while (pindex != vindex.end())
				{
					pa = a.begin() + *pindex;
					pd = d.begin() + *pindex;

					currentarray = *pa;
					currentdye = *pd;

					pa = a.begin();
					pd = d.begin();
					pc = c.begin();

					index = 0;

					while (pa!=a.end())
					{
						if (*pa == currentarray && *pd != currentdye)
						{
							pc += index;

							bool AddNew = true; // only add to tmpc vector when no this condition existing before.

							for (j=0;j<tmpc.size();j++)
							{
								if (tmpc.at(j) == *pc)
								{
									AddNew = false;
								}
							}

							if (AddNew)
							{
								tmpc.push_back(*pc);	
							}
						}

						index++;
						pa++;
						pd++;
					}
					pindex ++;
				}
			}

			// after the for loop, one block is found. 
			// store the block as an element in the designblock vector.
			// and delete the conditions involved from the unqiue condition vector, find next block.

			sort(tmpc.begin(),tmpc.end());
			vector <int>::iterator iter = unique(tmpc.begin(),tmpc.end());
			tmpc.erase(iter,tmpc.end());

			pa = a.begin();
			pc = c.begin();
			pd = d.begin();
			bool found = false;
			vector <int> blocka;
			vector <int> blockc;
			vector <int> blockd;
				
			// select the block according to the condition vector.

			while (pc!=c.end())
			{
				vector <int> tmp;  // vector tmp is just used in "includes" function. 
				tmp.push_back(*pc);
				found = includes(tmpc.begin(),tmpc.end(),tmp.begin(),tmp.end());

				if (found)
				{
					blocka.push_back(*pa);
					blockc.push_back(*pc);
					blockd.push_back(*pd);
				}

				pa++;
				pc++;
				pd++;					
			}

			// sort data according to array vector.

			vector <int> blocka_sorted(blocka.begin(),blocka.end());
			vector <int> blockc_sorted;
			vector <int> blockd_sorted;

			sort(blocka_sorted.begin(),blocka_sorted.end());

			vector <int>::iterator ps = blocka_sorted.begin();

			while (ps!=blocka_sorted.end())
			{
				for (i=0;i<blocka.size();i++)
				{
					if (blocka.at(i) == *ps)
					{
						blockc_sorted.push_back(blockc.at(i));
						blockd_sorted.push_back(blockd.at(i)); 
					}

				}

				ps += 2;				
			}			

			// make one block.

			Design block;
			block.setArray(blocka_sorted);
			block.setCond(blockc_sorted);
			block.setDye(blockd_sorted);

			designblock.push_back(block); 

			// erase the selected conditions, go for next block.

			pc_unique = c_unique.begin();
			vector <int> tmpc_unique;

			while (pc_unique != c_unique.end())
			{
				vector <int> tmp; // vector tmp is just used in "includes" function. 
				tmp.push_back(*pc_unique);

				found = includes(tmpc.begin(),tmpc.end(),tmp.begin(),tmp.end());
				
				if (!found)
				{
					tmpc_unique.push_back(*pc_unique);

				}

				pc_unique ++;
			}

			c_unique = tmpc_unique;
			pc_unique = c_unique.begin();
		}

	}else
	{
		designblock.push_back(des); 
	}

	return (designblock);
}
