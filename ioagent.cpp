#include "ioagent.h"

#include "scriptparser.h"
#include "vstring.h"
#include "pqmath.h"
#include "cleaner.h"

#include <dirent.h>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <iomanip>
#include <cmath>

using namespace std;

unsigned int *IOAgent::loadPQ(const char *srcfn, unsigned int &row, unsigned int&col)
{
    ifstream is(srcfn);
    if(!is.is_open())
    {
        cout<<"Fail to read "<<srcfn<<endl;
        exit(0);
    }

    is>>row>>col;
    unsigned int i, j;
    unsigned int *data = new unsigned[row*col];
    for(i = 0; i < row; i++)
    for(j = 0; j < col; j++)
    {
        is>>data[i*col+j];
    }
    is.close();
    return data;
}

float *IOAgent::loadPQInfo(const char *srcfn, unsigned int &pqn, unsigned int &pql, unsigned &pqm)
{
    ifstream is(srcfn);
    if(!is.is_open())
    {
        cout<<"Fail to read "<<srcfn<<endl;
        exit(0);
    }
    is>>pqn>>pql>>pqm;
    unsigned int i, j, k;
    float *data = new float[pqn*pql*pqm];

    for(i = 0; i < pqm; i++)
    {
        for(j = 0; j < pqn; j++)
        {
            for(k = 0; k < pql; k++)
            {
                is>>data[i*pqn*pql+j*pql+k];
            }
        }
    }
    is.close();
    return data;
}

float *IOAgent::loadMatrix(const char *srcfn, unsigned int &row, unsigned int &col)
{
    assert(srcfn);
    float vals[2] = {0};

    ifstream *inStrm = new ifstream();
    inStrm->open(srcfn, ios::in);
    if(inStrm->fail())
    {
        cout<<"Fail to read "<<srcfn<<endl;
        delete inStrm;
        exit(0);
    }

    (*inStrm)>>vals[0];
    (*inStrm)>>vals[1];
    row = (int)round(vals[0]);
    col = (int)round(vals[1]);
    float *mat = new float[row*col];
    memset(mat, 0, sizeof(float)*row*col);
    float *tmp = new float[col];
    unsigned int irow, idim, loc;

    for(irow = 0; irow < row; irow++)
    {
        loc = irow*col;
        for(idim = 0; idim < col; idim++)
        {
            (*inStrm) >>tmp[idim];
        }
        for(idim = 0; idim < col; idim++)
        {
            mat[loc + idim] = tmp[idim];
        }
        memset(tmp, 0, sizeof(float)*idim);
    }
    inStrm->close();
    assert(irow > 0);
    row = irow;
    delete inStrm;
    delete [] tmp;
    return mat;
}

float *IOAgent::load_fvecs(const char *fvecfn, unsigned int &d, unsigned int &r)
{

    float *mat = NULL, *vect = NULL, *ppmat = NULL;
    unsigned int di = 0, ir = 0;
    d = 0;
    r = 0;
    unsigned long bg = 0, fsize = 0, bfsz = 0;

    ifstream *inStrm = new ifstream(fvecfn, ios::in|ios::binary);
    if(!inStrm->is_open())
    {
        cout<<"File '"<<fvecfn<<"' cannot open for read!\n";
        return NULL;
    }

    //inStrm = new ifstream(fvecfn, ios::in|ios::binary);
    bg = (long)inStrm->tellg();
    inStrm->read((char*)&di,  sizeof(unsigned int));
    d = di;
    bfsz = d*sizeof(float);
    inStrm->seekg(0, ios::end);
    fsize = (unsigned long)inStrm->tellg() - bg;
    inStrm->close();
    r = fsize/(sizeof(unsigned int) + bfsz);
    if(r == 0)
    {
        cout<<"No data has been loaded!\n";
        r = d = 0;
        return NULL;
    }
    vect   = new float[d];
    mat    = new float[r*d];
    memset(mat, 0, sizeof(float)*r*d);
    ppmat  = mat;
    inStrm = new ifstream(fvecfn, ios::in|ios::binary);
    while(!inStrm->eof() && ir < r)
    {
        inStrm->read((char*)&di,  sizeof(int));

        if(di == 0)
            continue;

        bfsz = di*sizeof(float);
        inStrm->read((char*)vect, bfsz);
        memcpy(ppmat, vect, bfsz);
        memset(vect, 0, bfsz);
        ppmat = ppmat + d;
        di = 0;
        ir++;
    }

    delete [] vect;
    vect = NULL;
    inStrm->close();
    return mat;
}

float *IOAgent::loadDat(const char *fn, unsigned int &row, unsigned int &col)
{
    //TO-DO
    float *mat = NULL;
    string str;
    int id;
    float data;
    const char *loc = NULL;

    ifstream inStrm(fn);
    if(!inStrm.is_open())
    {
        cout<<"File '"<<fn<<"' cannot open for read!\n";
        exit(0);
    }

    getline(inStrm,str);
    loc = str.c_str();
    row = atoi(loc);
    loc = strstr(loc+1, " ");
    col = atoi(loc+1);
    mat = new float[row*col];
    memset(mat, 0, row*col*sizeof(float));

    for(unsigned i = 0; i < row; i++)
    {
        getline(inStrm,str);
        loc = str.c_str();
        int rowid = i*col;
        while(loc != NULL)
        {
            while(*loc == ' ')
                loc++;
            id = atoi(loc) ;
            loc = strstr(loc, " ");
            if(loc == NULL)
                break;
            while(*loc == ' ')
            {
                loc++;
            }
            data = atof(loc);
            loc = strstr(loc, " ");
            mat[rowid+id-1] = data;
        }
    }

    inStrm.close();
    return mat;
}

float *IOAgent::loadItms(const char *fn,    const char *idxKey,  unsigned int &row, unsigned int &col)
{
    ///  int ct = 0;
    ///  int datanum = 1000000;
    map<unsigned int, map<string, const char*>* >::iterator mit;
    map<unsigned int, map<string, const char*>* > itms;
    itms = ScriptParser::getItmMaps(fn);
    map<string, const char*> *crntItm;
    unsigned long  num = 0, dim = 0, i = 0, j = 0, idim = 0, loc = 0, irow = 0, vals[2];
    unsigned int ni = 0, di = 0, id = 0;
    ifstream inStrm;

    float *mat = NULL;
    float *tmp ;
    float data;
    string str;
    const char *ch  = NULL;
    bool Done = false;
    if(!strcmp(idxKey, "matrixdir"))
    {
        for(mit = itms.begin(); mit != itms.end(); mit++)
        {
            crntItm = mit->second;
            if(crntItm->find(idxKey) == crntItm->end())
            {
                continue;
            }
            vector<string> filenames;
            IOAgent::readFileName((*crntItm)[idxKey], filenames);
            mat = IOAgent::loadMatrix(filenames, row, col);
        }
        ScriptParser::clearItmMaps(itms);
        return mat;
    }
    else if(!strcmp(idxKey, "fvecsdir"))
    {
        for(mit = itms.begin(); mit != itms.end(); mit++)
        {
            crntItm = mit->second;
            if(crntItm->find(idxKey) == crntItm->end())
            {
                continue;
            }
            vector<string> filenames;
            IOAgent::readFileName((*crntItm)[idxKey], filenames);
            mat = IOAgent::load_fvecs(filenames, col, row);
        }
        ScriptParser::clearItmMaps(itms);
        return mat;
    }
    else if(!strcmp(idxKey, "matdir"))
    {
        for(mit = itms.begin(); mit != itms.end(); mit++)
        {
            crntItm = mit->second;
            if(crntItm->find(idxKey) == crntItm->end())
            {
                continue;
            }
            vector<string> filenames;
            IOAgent::readFileName((*crntItm)[idxKey], filenames);
            mat = IOAgent::loadDat(filenames, row, col);
        }
        ScriptParser::clearItmMaps(itms);
        return mat;
    }

    for(mit = itms.begin(); mit != itms.end(); mit++)
    {
        crntItm = mit->second;
        if(crntItm->find(idxKey) == crntItm->end())
        {
            cout<<"Parameter '"<<idxKey<<"' is missing!\n";
            exit(0);
        }
        Done = IOAgent::getMatrixInfo((*crntItm)[idxKey], ni, di);
        if(!Done)
        {
            cout<<"something goes wrong with '"<<(*crntItm)[idxKey]<<"'"<<endl;
        }
        ///cout<<ni<<"\t"<<di<<endl;
        assert(ni > 0);
        if(dim == 0)
        {
            dim = di;
        }
        assert(di == dim);
        num += ni;
    }
    row = num;
    col = dim;
    ///cout<<"row "<<row<<" col "<<col<<"\n";
    mat = new float[num*dim]; //if dealing with complex vector*2
    /// mat = new float[datanum*dim];
    tmp = new float[dim];
    loc = 0;
    for(mit = itms.begin(); mit != itms.end(); mit++)
    {
        crntItm = mit->second;
        if(crntItm->find(idxKey) == crntItm->end())
        {
            cout<<"Parameter '"<<idxKey<<"' is missing!\n";
            exit(0);
        }

        inStrm>>vals[0];
        inStrm>>vals[1];
        if(VString::endWith((*crntItm)[idxKey], ".txt"))
        {
            inStrm.open((*crntItm)[idxKey], ios::in);
            assert(inStrm.is_open());
            inStrm >>vals[0];
            irow = (int)round(vals[0]);

            for(j = 0; j < irow; j++)
            {

                for(idim = 0; idim < dim; idim++)
                {
                    inStrm >>tmp[idim];
                }
                for(idim = 0; idim < dim; idim++)
                {
                    mat[loc + idim] = tmp[idim];
                }
                memset(tmp, 0, sizeof(float)*idim);
                loc += dim;
                /**     ct++;
                  if(ct == datanum)
                     {
                         delete [] tmp;
                         tmp = NULL;
                         ScriptParser::clearItmMaps(itms);
                         row = datanum;
                         return mat;

                     }
                 /**/
            }
        }
        else if(VString::endWith((*crntItm)[idxKey], ".fvecs"))
        {
            inStrm.open((*crntItm)[idxKey], ios::in|ios::binary);

            assert(inStrm.is_open());
            unsigned long bg = inStrm.tellg();
            inStrm.read((char*)&dim, sizeof(int));
            inStrm.seekg(0, ios::end);
            unsigned long sz = inStrm.tellg();
            num = (sz - bg)/(sizeof(int) + col*sizeof(float));
            inStrm.seekg(0, ios::beg);
            j = 0;
            while(!inStrm.eof() && j < num)
            {
                inStrm.read((char*)&dim,  sizeof(int));
                if(dim == 0)
                    continue;

                inStrm.read((char*)tmp, dim*sizeof(float));

                memcpy(mat+loc, tmp, dim*sizeof(float));
                memset(tmp, 0, dim*sizeof(float));
                loc += dim;
                dim = 0;
                j++;
                /**    ct++;
                    if(ct == datanum)
                    {
                        delete [] tmp;
                        tmp = NULL;
                        ScriptParser::clearItmMaps(itms);
                        row = datanum;
                        return mat;

                    }
                 /**/
            }
        }
        else if(VString::endWith((*crntItm)[idxKey], ".mat"))
        {
            inStrm.open((*crntItm)[idxKey], ios::in);
            assert(inStrm.is_open());
            getline(inStrm,str);
            ch = str.c_str();
            irow += atoi(ch);

            for(i = 0; i < irow; i++)
            {
                getline(inStrm,str);
                ch = str.c_str();
                while(ch != NULL)
                {
                    while(*ch == ' ')
                        ch++;
                    id = atoi(ch) ;
                    ch = strstr(ch, " ");
                    if(ch == NULL)
                    break;
                    while(*ch == ' ')
                        ch++;
                    data = atof(ch);
                    ch = strstr(ch, " ");
                    mat[loc+id-1] = data;
                }
                loc += col;
                /**     ct++;
                     if(ct == datanum)
                     {
                         delete [] tmp;
                         tmp = NULL;
                         ScriptParser::clearItmMaps(itms);
                         row = datanum;
                         return mat;

                     }
                 /**/
            }
        }
        ///filling your code
        inStrm.close();
    }

    delete [] tmp;
    tmp = NULL;
    ScriptParser::clearItmMaps(itms);

    return mat;
}


float *IOAgent::loadItms(const char *fn,    const char *idxKey, const unsigned int line,  unsigned int &row, unsigned int &col)
{
    map<unsigned int, map<string, const char*>* >::iterator mit;
    map<unsigned int, map<string, const char*>* > itms;
    itms = ScriptParser::getItmMaps(fn);
    map<string, const char*> *crntItm;
    unsigned long  num = 0, dim = 0, i = 0, j = 0, idim = 0, loc = 0, irow = 0, vals[2];
    unsigned int ni = 0, di = 0, id = 0;
    ifstream inStrm;

    float *mat = NULL;
    float *tmp ;
    float data;
    string str;
    const char *ch  = NULL;
    bool Done = false;
    if(!strcmp(idxKey, "matrixdir"))
    {
        for(mit = itms.begin(); mit != itms.end(); mit++)
        {
            crntItm = mit->second;
            if(crntItm->find(idxKey) == crntItm->end())
            {
                continue;
            }
            vector<string> filenames;
            IOAgent::readFileName((*crntItm)[idxKey], filenames);
            mat = IOAgent::loadMatrix(filenames, row, col);
        }
        ScriptParser::clearItmMaps(itms);
        return mat;
    }
    else if(!strcmp(idxKey, "fvecsdir"))
    {
        for(mit = itms.begin(); mit != itms.end(); mit++)
        {
            crntItm = mit->second;
            if(crntItm->find(idxKey) == crntItm->end())
            {
                continue;
            }
            vector<string> filenames;
            IOAgent::readFileName((*crntItm)[idxKey], filenames);
            mat = IOAgent::load_fvecs(filenames, col, row);
        }
        ScriptParser::clearItmMaps(itms);
        return mat;
    }
    else if(!strcmp(idxKey, "matdir"))
    {
        for(mit = itms.begin(); mit != itms.end(); mit++)
        {
            crntItm = mit->second;
            if(crntItm->find(idxKey) == crntItm->end())
            {
                continue;
            }
            vector<string> filenames;
            IOAgent::readFileName((*crntItm)[idxKey], filenames);
            mat = IOAgent::loadDat(filenames, row, col);
        }
        ScriptParser::clearItmMaps(itms);
        return mat;
    }

    for(mit = itms.begin(); mit != itms.end(); mit++)
    {
        crntItm = mit->second;
        if(crntItm->find(idxKey) == crntItm->end())
        {
            cout<<"Parameter '"<<idxKey<<"' is missing!\n";
            exit(0);
        }
        Done = IOAgent::getMatrixInfo((*crntItm)[idxKey], ni, di);
        if(!Done)
        {
            cout<<"something goes wrong with '"<<(*crntItm)[idxKey]<<"'"<<endl;
        }
        ///cout<<ni<<"\t"<<di<<endl;
        assert(ni > 0);
        if(dim == 0)
        {
            dim = di;
        }
        assert(di == dim);
        num += ni;
    }
    row = num;
    col = dim;

    if(row < line)
    {
        cout<<"erro not enough line row = "<<row<<" \n";
        exit(0);
    }

    ///cout<<"row "<<row<<" col "<<col<<"\n";
    mat = new float[line*dim]; //if dealing with complex vector*2
    tmp = new float[dim];
    loc = 0;
    row = 0;
    for(mit = itms.begin(); mit != itms.end(); mit++)
    {
        crntItm = mit->second;
        if(crntItm->find(idxKey) == crntItm->end())
        {
            cout<<"Parameter '"<<idxKey<<"' is missing!\n";
            exit(0);
        }

        inStrm>>vals[0];
        inStrm>>vals[1];
        if(VString::endWith((*crntItm)[idxKey], ".txt"))
        {
            inStrm.open((*crntItm)[idxKey], ios::in);
            assert(inStrm.is_open());
            inStrm >>vals[0];
            irow = (int)round(vals[0]);

            for(j = 0; j < irow; j++)
            {

                for(idim = 0; idim < dim; idim++)
                {
                    inStrm >>tmp[idim];
                }
                for(idim = 0; idim < dim; idim++)
                {
                    mat[loc + idim] = tmp[idim];
                }
                memset(tmp, 0, sizeof(float)*idim);
                loc += dim;
                row ++;
                if(row == line)
                {
                    inStrm.close();
                    delete [] tmp;
                    tmp = NULL;
                    ScriptParser::clearItmMaps(itms);
                    return mat;
                }
            }
        }
        else if(VString::endWith((*crntItm)[idxKey], ".fvecs"))
        {
            inStrm.open((*crntItm)[idxKey], ios::in|ios::binary);

            assert(inStrm.is_open());
            unsigned long bg = inStrm.tellg();
            inStrm.read((char*)&dim, sizeof(int));
            inStrm.seekg(0, ios::end);
            unsigned long sz = inStrm.tellg();
            num = (sz - bg)/(sizeof(int) + col*sizeof(float));
            inStrm.seekg(0, ios::beg);
            j = 0;
            while(!inStrm.eof() && j < num)
            {
                inStrm.read((char*)&dim,  sizeof(int));
                if(dim == 0)
                    continue;

                inStrm.read((char*)tmp, dim*sizeof(float));

                memcpy(mat+loc, tmp, dim*sizeof(float));
                memset(tmp, 0, dim*sizeof(float));
                loc += dim;
                dim = 0;
                j++;
                row ++;
                if(row == line)
                {
                    inStrm.close();
                    delete [] tmp;
                    tmp = NULL;
                    ScriptParser::clearItmMaps(itms);
                    return mat;
                }
            }
        }
        else if(VString::endWith((*crntItm)[idxKey], ".mat"))
        {
            inStrm.open((*crntItm)[idxKey], ios::in);
            assert(inStrm.is_open());
            getline(inStrm,str);
            ch = str.c_str();
            irow += atoi(ch);

            for(i = 0; i < irow; i++)
            {
                getline(inStrm,str);
                ch = str.c_str();
                while(ch != NULL)
                {
                    while(*ch == ' ')
                        ch++;
                    id = atoi(ch) ;
                    ch = strstr(ch, " ");
                    if(ch == NULL)
                        break;
                    while(*ch == ' ')
                        ch++;
                    data = atof(ch);
                    ch = strstr(ch, " ");
                    mat[loc+id-1] = data;
                }
                loc += col;
                row ++;
                if(row == line)
                {
                    inStrm.close();
                    delete [] tmp;
                    tmp = NULL;
                    ScriptParser::clearItmMaps(itms);
                    return mat;
                }
            }
        }
        ///filling your code
        inStrm.close();
    }

    delete [] tmp;
    tmp = NULL;
    ScriptParser::clearItmMaps(itms);

    return mat;
}

SparseMatrix IOAgent::loadBOVW(const char * fn, unsigned int &row, unsigned int &col)
{
    SparseMatrix sdata;
    fstream is(fn);
    ///  char sparsefn[1024];
    /// sprintf(sparsefn, "s%s", fn);
    ///fstream os(sparsefn);
    if(!is.is_open())
    {
        cout<<"file '"<<fn<<"' can not open for read\n";
        exit(0);
    }

    int n, tmp1, tmp2, nozore;

    string str;
    n = 0;
    nozore = 0;
    while(is>>tmp1)
    {
        nozore += tmp1;
        n++;
        getline(is, str);
    }

    row = n;

    double *norm;

    norm = new double[row];
    memset(norm, 0, sizeof(double));
    sdata.data = new float[nozore];
    sdata.index = new unsigned int[nozore];
    sdata.col = new unsigned int[row + 1];
    memset(sdata.data, 0, nozore*sizeof(float));
    memset(sdata.index, 0, nozore*sizeof(unsigned int));
    memset(sdata.col, 0, (row+1)*sizeof(unsigned int));

    is.clear();
    is.seekg(0, ios::beg);

    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;
    col = 0;
    while(is>>n)
    {
        sdata.col[i] = k;
        for(j = 0; j < n; j++)
        {
            is>>tmp1>>tmp2;
            sdata.data[k] = tmp2;
            sdata.index[k] = tmp1;
            if(tmp1 > col)
                col = tmp1;
            k++;
            norm[i] += tmp2*tmp2;
        }
        i++;
    }
    col++;
    sdata.col[i] = k;

    /**norm**/

    for(i = 0; i < row; i++)
    {
        for(j = sdata.col[i]; j < sdata.col[i+1]; j++)
        {
            sdata.data[j] /= norm[i];
        }
    }


    is.close();
    delete [] norm;
    return sdata;
}

unsigned int IOAgent::getMatrixInfo(const char *matfn, unsigned int &row, unsigned int &col)
{
    ifstream inStrm;
    row = 0;
    const char *loc = NULL;
    string str;
    int tmpDim = 0;
    if(VString::endWith(matfn, ".fvecs"))
    {
        inStrm.open(matfn, ios::in|ios::binary);
        unsigned long bg = inStrm.tellg();
        inStrm.read((char*)&tmpDim, sizeof(int));
        inStrm.seekg(0, ios::end);
        unsigned long sz = inStrm.tellg();
        col = tmpDim;
        row = (sz - bg)/(sizeof(int) + col*sizeof(float));
    }
    else if(VString::endWith(matfn, ".txt"))
    {
        inStrm.open(matfn, ios::in);
        inStrm>>row>>col;
    }
    else if(VString::endWith(matfn, ".dat"))
    {
        inStrm.open(matfn, ios::in);
        if(!inStrm.is_open())
        {
            cout<<"File "<<matfn<<" cannot open for read!\n";
            exit(0);
        }

        getline(inStrm,str);
        loc = str.c_str();
        row = atoi(loc);
        loc = strstr(loc+1, " ");
        col = atoi(loc+1);
    }
    inStrm.close();
    return row;
}

float *IOAgent::loadMatrix(vector<string> &filenames, unsigned int &row, unsigned int &col)
{
    float vals[2] = {0};
    int i = 0, j = 0;
    row = 0;
    col = 0;

    ifstream *inStr = new ifstream();

    for(i = 0; i < filenames.size(); i++)
    {
        inStr->open(filenames[i], ios::in);
        if(!inStr->is_open())
        {
            cout<<"Fail to read "<<filenames[i]<<endl;
            delete inStr;
            exit(0);
        }
        (*inStr)>>vals[0];
        (*inStr)>>vals[1];
        col = (int)round(vals[1]);
        row += (int)round(vals[0]);
        inStr->close();
    }
    float *mat = new float[row*col];
    memset(mat, 0, sizeof(float)*row*col);
    float *tmp = new float[col];
    unsigned int irow, idim, loc;

    loc = 0;
    for(i = 0; i < filenames.size(); i++)
    {

        inStr->open(filenames[i], ios::in);
        if(inStr->fail())
        {
            cout<<"Fail to read "<<filenames[i]<<endl;
            delete inStr;
            exit(0);
        }

        (*inStr)>>vals[0];
        (*inStr)>>vals[1];

        irow = (int)round(vals[0]);

        for(j = 0; j < irow; j++)
        {
            for(idim = 0; idim < col; idim++)
            {
                (*inStr) >>tmp[idim];
            }
            for(idim = 0; idim < col; idim++)
            {
                mat[loc + idim] = tmp[idim];
            }
            memset(tmp, 0, sizeof(float)*idim);
            loc += col;
        }
        inStr->close();
        assert(irow > 0);
    }

    delete inStr;
    delete [] tmp;
    return mat;
}

float *IOAgent::load_fvecs(vector<string> &filenames, unsigned int &d, unsigned int &r)
{

    float *mat = NULL, *vect = NULL, *ppmat = NULL;
    unsigned int di = 0, ir = 0, i = 0;
    d = 0;
    unsigned long bg = 0, fsize = 0, bfsz = 0;
    unsigned int *tmpr = NULL;
    tmpr = new unsigned int[filenames.size()];
    ifstream *inStrm = new ifstream();

    r = 0;
    for(i = 0; i <  filenames.size(); i++)
    {

        inStrm->open(filenames[i], ios::in|ios::binary);
        if(!inStrm->is_open())
        {
            cout<<"File '"<<filenames[i]<<"' cannot open for read!\n";
            return NULL;
        }
        bg = (long)inStrm->tellg();
        inStrm->read((char*)&di,  sizeof(unsigned int));
        d = di;
        bfsz = d*sizeof(float);
        inStrm->seekg(0, ios::end);
        fsize = (unsigned long)inStrm->tellg() - bg;
        inStrm->close();
        tmpr[i] = fsize/(sizeof(unsigned int) + bfsz);
        r += tmpr[i];
    }
    if(r == 0)
    {
        cout<<"No data has been loaded!\n";
        r = d = 0;
        return NULL;
    }
    vect   = new float[d];
    mat    = new float[r*d];
    memset(mat, 0, sizeof(float)*r*d);
    ppmat  = mat;
    for(i = 0; i < filenames.size(); i++)
    {
        inStrm->open(filenames[i], ios::in|ios::binary);
        ir = 0;
        while(!inStrm->eof() && ir < tmpr[i])
        {
            inStrm->read((char*)&di,  sizeof(int));

            if(di == 0)
                continue;

            bfsz = di*sizeof(float);
            inStrm->read((char*)vect, bfsz);
            memcpy(ppmat, vect, bfsz);
            memset(vect, 0, bfsz);
            ppmat = ppmat + d;
            di = 0;
            ir++;
        }
        inStrm->close();
    }
    delete [] vect;
    vect = NULL;
    return mat;
}

float *IOAgent::loadDat(vector<string> &filenames, unsigned int &row, unsigned int &col)
{
    //TO-DO
    float *mat = NULL;
    string str;
    unsigned int id, j, i, rowid, irow;
    float data;
    const char *loc = NULL;

    ifstream inStrm;
    for(i = 0; i < filenames.size(); i++)
    {
        inStrm.open(filenames[i], ios::in);
        if(!inStrm.is_open())
        {
            cout<<"File "<<filenames[i]<<" cannot open for read!\n";
            exit(0);
        }

        getline(inStrm,str);
        loc = str.c_str();
        row += atoi(loc);
        loc = strstr(loc+1, " ");
        col = atoi(loc+1);
        inStrm.close();
    }
    mat = new float[row*col];
    memset(mat, 0, row*col*sizeof(float));

    rowid = 0;
    for(j = 0; j < filenames.size(); j++)
    {
        inStrm.open(filenames[i], ios::in);
        getline(inStrm,str);
        loc = str.c_str();
        irow += atoi(loc);
        {
            getline(inStrm,str);
            loc = str.c_str();
            while(loc != NULL)
            {
                while(*loc == ' ')
                    loc++;
                id = atoi(loc) ;
                loc = strstr(loc, " ");
                while(*loc == ' ')
                    loc++;
                data = atof(loc);
                loc = strstr(loc, " ");
                mat[rowid+id-1] = data;
            }
            rowid += col;
        }
        inStrm.close();
    }
    inStrm.close();
    return mat;
}


void IOAgent::saveDat(const char* fn, const unsigned int& row, const unsigned int &col, float *Dat)
{
    unsigned int i, j;
    ofstream ofStrm(fn);
    if(!ofStrm.is_open())
    {
        cout<<"File "<<fn<<" cannot open for write!\n";
        exit(0);
    }

    int counter = 0;
    for(i = 0; i < row; i++)
    {
        for(j = 0; j < col; j++)
        {
            if(Dat[i*col + j] != 0)
                counter++;
        }
    }

    ofStrm<<row<<" "<<col<<" "<<counter<<"\n";

    for(i = 0; i < row; i++)
    {
        for(j = 0; j < col; j++)
        {
            if(Dat[i*col + j] != 0)
            {
                ofStrm<<" "<<(j+1)<<" "<<Dat[i*col + j];
            }
        }
        ofStrm<<"\n";
    }

    ofStrm.close();
}

bool IOAgent::readFileName(const char *srcdir,vector< string > &filename)
{
    DIR *srcdp = NULL;
    struct dirent *result = NULL, entry;
    if((srcdp  = opendir(srcdir)) == NULL)
    {
        cout<<"Source directory '"<<srcdir<<"' do not exist!\n";
        return false;
    }
    while (readdir_r(srcdp, &entry, &result) == 0 && result != NULL)
    {
        //check each entry.d_name by yourself
        if(!strcmp(entry.d_name, ".") || !strcmp(entry.d_name, ".."))
            continue;
        filename.push_back(entry.d_name);

    }
    return true;
}

float *IOAgent::loadCSV(const char *srcFn, const char delm, unsigned int &row, unsigned int &col)
{
     float *mat = NULL;
     unsigned int count = 0, idx = 0, i = 0;
     vector<float> vals;
     row = 0;  col = 0;
     string txt;
     ifstream *inStrm = new ifstream(srcFn, ios::in);
     while(inStrm->good())
     {
         getline(*inStrm, txt, '\n');
         if(col == 0)
         {
            col = VString::countsof(txt, delm) + 1;
         }
         if(txt.length() > 1)
            count++;
     }
     ///cout<<"col: "<<col<<endl;
     if(count == 0)
     {
         return NULL;
     }

     row = count;
     mat = new float[row*col];
     inStrm->clear();
     inStrm->seekg(0, ios::beg);
     while(inStrm->good())
     {
         getline(*inStrm, txt, '\n');
         if(txt.size() > 1)
         {
            VString::str2float(txt.c_str(), ',', vals);
            for(i = 0; i < vals.size(); i++, idx++)
            {
               mat[idx] = vals[i];
            }
            vals.clear();
            txt.clear();
         }
      }
      ///cout<<mat[idx-1]<<endl;

      inStrm->close();
      return mat;
}

void IOAgent::saveEx(vector<string> &filenames, float *e, float *p, const char *exfile)
{
    ofstream oStr(exfile);
    if(!oStr.is_open())
    {
        cout<<"\ncant not open file "<<exfile<<" for write"<<endl;
        exit(0);
    }
    oStr<<setw(8)<<"dataset"<<"\t"<<setw(8)<<"k=5"<<"\t"<<setw(8)<<"k=10"<<"\t"<<setw(8)<<"k=15"<<"\t"<<setw(8)<<"k=20"<<"\t";
    oStr<<setw(8)<<"p=5"<<"\t"<<setw(8)<<"p=10"<<"\t"<<setw(8)<<"p=15"<<"\t"<<setw(8)<<"p=20"<<"\t"<<endl;

    int s =filenames.size();
    int id = 0;
    for(int i = 0; i < s; i++)
    {
        oStr<<setw(8)<<filenames[i]<<"\t"<<setw(8)<<e[id+0]<<"\t"<<setw(8)<<e[id+1]<<"\t"<<setw(8)<<e[id+2]<<"\t"<<setw(8)<<e[id+3]<<"\t";
        oStr<<setw(8)<<p[id+0]<<"\t"<<setw(8)<<p[id+1]<<"\t"<<setw(8)<<p[id+2]<<"\t"<<setw(8)<<p[id+3]<<endl;
        id += 4;

    }
    oStr.close();
}

SparseMatrix IOAgent::loadSparse(const char *fn, unsigned int &row, unsigned int &col)
{
    SparseMatrix sdata;
    string str;
    unsigned int id, i, nozore, j;
    float data;
    const char *loc = NULL;

    ifstream inStrm(fn);
    if(!inStrm.is_open())
    {
        cout<<"File '"<<fn<<"' cannot open for read!\n";
        exit(0);
    }

    getline(inStrm,str);
    loc = str.c_str();
    row = atoi(loc);
    loc = strstr(loc+1, " ");
    col = atoi(loc+1);
    loc = strstr(loc+1, " ");
    nozore = atoi(loc+1);
    sdata.data = new float[nozore];
    sdata.index = new unsigned int[nozore];
    sdata.col = new unsigned int[row + 1];
    memset(sdata.data, 0, nozore*sizeof(float));
    memset(sdata.index, 0, nozore*sizeof(unsigned int));
    memset(sdata.col, 0, (row+1)*sizeof(unsigned int));

    cout<<"row\t"<<row<<"\tcol\t"<<col<<"\n";

    j = 0;
    for(i = 0; i < row; i++)
    {
        getline(inStrm,str);
        loc = str.c_str();
        sdata.col[i] = j;
        while(loc != NULL)
        {
            while(*loc == ' ')
                loc++;
            id = atoi(loc) ;
            loc = strstr(loc, " ");
            while(*loc == ' ')
                loc++;
            data = atof(loc);
            loc = strstr(loc, " ");
            sdata.data[j] = data;
            sdata.index[j] = id - 1;
            j++;
        }

        if(i == row -1)
            sdata.col[row] = j;
    }

    inStrm.close();
    return sdata;
}

void IOAgent::addRslt(const char *filename, float *e, float *p, const char *result)
{
    ofstream oStr(result, ios::out | ios::app);
    if(!oStr.is_open())
    {
        cout<<"\ncant not open file "<<result<<" for write"<<endl;
        exit(0);
    }
    oStr<<setw(8)<<filename<<"\t"<<setw(8)<<e[0]<<"\t"<<setw(8)<<e[1]<<"\t"<<setw(8)<<e[2]<<"\t"<<setw(8)<<e[3]<<"\t";
    oStr<<setw(8)<<p[0]<<"\t"<<setw(8)<<p[1]<<"\t"<<setw(8)<<p[2]<<"\t"<<setw(8)<<p[3]<<endl;
    oStr.close();
}

void IOAgent::saveAsFullMat(const float* dat, unsigned int &row, unsigned int &col, const char *file)
{
    ofstream oStr(file);
    if(!oStr.is_open())
    {
        cout<<"file "<<file<<" can not open for write full mat\n";
        exit(0);
    }
    unsigned int loc = 0, i = 0, j = 0;
    for(i = 0; i < row; i++)
    {
        for(j = 0; j < col; j++)
        {
            oStr<<dat[loc + j]<<" ";
        }
        loc += col;
        oStr<<"\n";
    }
    oStr.close();
}

void IOAgent::saveSparse(const char* fn, const unsigned int& row, const unsigned int &col, const SparseMatrix &sdata)
{
    ofstream ofStrm(fn);
    if(!ofStrm.is_open())
    {
        cout<<"File '"<<fn<<"' cannot open for write!\n";
        exit(0);
    }

    int counter = sdata.col[row];


    ofStrm<<row<<" "<<col<<" "<<counter<<"\n";

    for(unsigned i = 0; i < row; i++)
    {
        for(unsigned j = sdata.col[i]; j < sdata.col[i+1]; j++)
        {
            ofStrm<<" "<<sdata.index[j]+1<<" "<<sdata.data[j];
        }
        ofStrm<<"\n";
    }

    ofStrm.close();
}

void IOAgent::loadClust(const char *fn, map<unsigned int, set<int> *> &cluster)
{
/// colnum1 data11 data12 .....datacolnum1
/// colnum2 data21 data22....
    ifstream is(fn);
    if(!is.is_open())
    {
        cout<<"file "<<fn<<" can not opened for read"<<"\n";
        exit(0);
    }

    int tmp, tmp1, i, j;
    set<int> *tmpset;

    j = 0;
    while(is>>tmp)
    {
        j++;
        tmpset = new set<int>;
        for(i = 0; i < tmp; i++)
        {
            is>>tmp1;
            tmpset->insert(tmp1);
        }

        cluster.insert(pair<int, set<int>* >(j, tmpset));

    }
    is.close();
}

void IOAgent::loadmClust(const char *srcFn, map<unsigned int, set<int> *> &clusts)
{
///label to clust
///label1
///label2
///.....
///labeln

    string tmp;
    vector<string> clust;
    set<string> temp;
    set<string>::iterator it;
    map<string, unsigned int> mp;
    int i;

    if(clusts.size())
    {
        cout<<"clear clusts *******"<<endl;
        Cleaner::clearClust(clusts);
    }

    ifstream inStr(srcFn);

    if(!inStr.is_open())
    {
        cout<<"Evaluator::loadClust open file "<<srcFn<<" fail"<<endl;
        exit(0);
    }
    while(inStr>>tmp)
    {
        clust.push_back(tmp);
    }

    inStr.close();


    int s = clust.size();
    for(int i = 0; i < s; i++)
    {
        temp.insert(clust[i]);
    }

    i = 0;
    for(it = temp.begin(); it != temp.end(); it++)
    {
        set<int > * tmpclust = new set<int>;
        mp.insert(pair<string,unsigned int>(*it, i));
        clusts.insert(pair<unsigned int, set<int> *>(i, tmpclust));
        i++;

    }

    for(int i = 0; i < s; i++)
    {
        clusts[mp[clust[i]]]->insert(i);
    }

    clust.clear();

    return ;
}

void IOAgent::saveClust(const char *fn, map<unsigned int, set<int> *> &cluster)
{
    ofstream os(fn);
    if(!os.is_open())
    {
        cout<<"file "<<fn<<" can not opened for write"<<"\n";
        exit(0);
    }

    int i, n, s;
    set<int>::iterator it;
    set<int> *pset;
    n = cluster.size();
    for(i = 0; i < n; i++)
    {
        pset = cluster[i];
        s = pset->size();
        os<<s<<"\t";
        for(it = pset->begin(); it != pset->end(); it++)
        {
            os<<*it<<"\t";
        }
        os<<"\n";
    }

    os.close();
}

void IOAgent::test()
{
    unsigned int row = 0, col = 0;
    float *data = NULL;

    const char *srcfn = "/home/wlzhao/datasets/clust/SUSY.csv";
    const char *dstfn = "/home/wlzhao/datasets/clust/SUSY_5m.txt";
    data = IOAgent::loadCSV(srcfn, ',', row, col);
    cout<<row<<"\t"<<col<<endl;
    delete [] data;

}
