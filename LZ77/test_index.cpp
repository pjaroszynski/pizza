#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "static_selfindex.h"

int main(int argc,char** argv){
    std::cout << "Creating Index"<<std::endl;
    lz77index::static_selfindex* idx = lz77index::static_selfindex_lz77::build(argv[1]);
    if(idx==NULL){
        std::cerr<<"Index is NULL"<<std::endl;
        return EXIT_FAILURE;
    }
    unsigned int count = 0;
    unsigned char* pattern = new unsigned char[strlen(argv[2])+1];
    strcpy((char*)pattern,argv[2]);
    //unsigned char pattern[] = {'e','x','p','e','r','i','m','e','n','t','\0'};
    //unsigned char pattern[] = {'U','l','m','\0'};
    //unsigned char pattern[] = {'a','\0'};
    std::vector<unsigned int>* pos = idx->locate(pattern,&count);
    std::cout<<"Found: "<<count<<" occurrences"<<std::endl;
    for(unsigned int i=0;i<count;i++){
        std::cout<<(*pos)[i]<<" ";
    }
    std::cout<<std::endl;
    if(pos!=NULL){
        delete [] pos;
    }
    delete [] pattern;/*
    for(unsigned int i=0;i<20;i++){
        unsigned char* substring = idx->display(i,i);
        std::cout<<substring;
        delete [] substring;
    }
    std::cout<<std::endl;*/
    //unsigned char* substring = idx->display(0,2000000);
    //std::cout<<substring<<std::endl;
    //delete [] substring;
    delete idx;
    return EXIT_SUCCESS;
}
/*int main2(int argc,char** argv){
    FILE* fp;    
    unsigned int count=0;
    unsigned int len=strlen(argv[1]);
    unsigned char* patt1 = new unsigned char[len+1];
    strncpy((char*)patt1,argv[2],len+1);
    unsigned char* patt2 = new unsigned char[len+1];
    strncpy((char*)patt2,argv[2],len+1);
    unsigned char* patt3 = new unsigned char[len+1];
    strncpy((char*)patt3,argv[2],len+1);
    unsigned char* patt4 = new unsigned char[len+1];
    strncpy((char*)patt4,argv[2],len+1);
    static_selfindex_none* idx;
    std::cout << "Creating Index"<<std::endl;
    idx = new static_selfindex_none((unsigned char*)argv[1]);
    std::cout<<"Size: "<<idx->size()<<std::endl;    
    std::cout<<"Finding Occurrences"<<std::endl;
    idx->locate(patt1,&count);
    std::cout<<std::endl;    
    std::cout<<std::endl<<count<<" "<<idx->count(patt2)<<std::endl;
    unsigned char* str;
    for(unsigned int i=0;i<len;i++){
        str=idx->display(i,i);
        std::cout<<str<<std::flush;
        delete [] str;
    }
    std::cout<<std::endl;
    str=idx->display(0,len-1);
    std::cout<<str<<std::endl;
    delete [] str;    
    fp = fopen("text.none.idx","w");
    unsigned int ret = idx->save(fp);
    std::cout<<"Save:"<<ret<<std::endl;
    fclose(fp);
    delete idx;
    std::cout << "Loading Index"<<std::endl;
    fp = fopen("text.none.idx","r");
    idx = static_selfindex_none::load(fp);
    fclose(fp);
    std::cout<<"Size: "<<idx->size()<<std::endl;
    std::cout << "Finding Occurrences"<<std::endl;
    idx->locate(patt3,&count);
    std::cout<<std::endl;    
    std::cout<<std::endl<<count<<" "<<idx->count(patt4)<<std::endl;
    for(unsigned int i=0;i<len;i++){
        str=idx->display(i,i);
        std::cout<<str<<std::flush;
        delete [] str;
    }
    std::cout<<std::endl;
    str=idx->display(0,len-1);
    std::cout<<str<<std::endl;
    delete [] str;    
    delete idx;
    delete [] patt1; delete [] patt2; delete [] patt3; delete [] patt4;
    return EXIT_SUCCESS;
}*/

