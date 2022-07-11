#ifndef NNITEM_H
#define NNITEM_H

#include <algorithm>
#include <vector>
#include <list>
#include <map>

using namespace std;

struct PairItm
{
public :
    unsigned int idx;
    int numb;
    double val;
    static int LLcomparer(const PairItm &a, const PairItm &b)
    {
        if(a.numb < b.numb)
        {
            return 1;
        }else if(a.numb == b.numb)
        {
            if(a.val < b.val)
            {
               return 1;
            }else{
               return 0;
            }
        }else{
           return 0;
        }
    }

    static int LGcomparer(const PairItm &a, const PairItm &b)
    {
        if(a.numb > b.numb)
        {
            return 1;
        }else if(a.numb == b.numb)
        {
            if(a.val < b.val)
            {
               return 1;
            }else{
               return 0;
            }
        }else{
           return 0;
        }
    }
};

struct MiniNN
{
public:
    unsigned int idx;
    float val;
    unsigned char nw;

    /**ascending order**/
    static int LLcomparer(const MiniNN &a, const MiniNN &b)
    {
        return (a.val < b.val);
    }



    /**descending order**/
    static int LGcomparer(const MiniNN &a, const MiniNN &b)
    {
        return (a.val > b.val);
    }


};

struct NbHood{
   vector<unsigned int> oldnb;
   vector<unsigned int> newnb;
   vector<unsigned int> roldnb;
   vector<unsigned int> rnewnb;
   unsigned short rnew, rold;
};

class NNItem
{
public:
    NNItem(const unsigned int index, const double dist)
    {
        this->index = index;
        this->val   = dist;
        this->size  = 0;
        this->dvd   = true;
    }
    NNItem(const unsigned int index, const double dist, const int sz)
    {
        this->index = index;
        this->val   = dist;
        this->size  = sz;
        this->dvd   = true;
    }
    unsigned int index;
    double val;
    int   size;
    bool  dvd;

    /**ascending order**/
    static int LLcomparer(const NNItem *a,const NNItem *b)
    {
        return (a->val < b->val);
    }

    /**descending order**/
    static int LGcomparer(const NNItem *a, const NNItem *b)
    {
        return (a->val > b->val);
    }

    /**ascending order**/
    static int LLIDXcomparer(const NNItem *a, const NNItem *b)
    {

        return (a->index < b->index);
    }

    /**ascending order, consider the priority first**/
    static int LLVALcomparer(const NNItem *a, const NNItem *b)
    {
        if((a->dvd == false) && (b->dvd == false))
        {
            return false;
        }
        else if((a->dvd == false) && b->dvd)
        {
            return false;
        }
        else if(a->dvd && (b->dvd == false))
        {
            return true;
        }
        else
        {
            if(a->val < b->val)
            {
                return true;
            }
            else if(a->val == b->val)
            {
                if(a->size > b->size)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }

    /**descending order, consider the priority first**/
    static int LGVALcomparer(const NNItem *a, const NNItem *b)
    {
        if((a->dvd == false) && (b->dvd == false))
        {
            return true;
        }
        else if((a->dvd == false) && b->dvd)
        {
            return false;
        }
        else if(a->dvd && (b->dvd == false))
        {
            return true;
        }
        else
        {
            if(a->val > b->val)
            {
                return true;
            }
            else if(a->val == b->val)
            {
                if(a->size > b->size)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }

    /**descending order, consider the priority first**/
    static int LGSZcomparer(const NNItem *a, const NNItem *b)
    {
        if((a->dvd == false) && (b->dvd == false))
        {
            return true;
        }
        else if((a->dvd == false) && b->dvd)
        {
            return false;
        }
        else if(a->dvd && (b->dvd == false))
        {
            return true;
        }
        else
        {
            if(a->size > b->size)
            {
                return true;
            }
            else if(a->size == b->size)
            {
                if(a->val > b->val)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }
};


/******************************************/

#endif
