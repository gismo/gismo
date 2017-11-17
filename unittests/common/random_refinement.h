
#pragma once




//lvl- max level of domain refinement, max number of boxes inserted
//per level, basis, alignement with the prewious level
std::vector<unsigned int> random_refinement(int lvl, int max_nb, 
                                            gismo::gsHTensorBasis<2> *bas, 
                                            bool aligned = false)
{
    using namespace gismo;

    //int max_x, max_y;
    int boxes_in_lvl, bi,span_size;
    std::vector<unsigned int> q;
    std::vector<gsVector<unsigned int> > boxes;// boxes stores the boxes inserted to the previous level
    std::vector<gsVector<unsigned int> > boxes_new;
    gsVector<unsigned int, 2> i1;
    gsVector<unsigned int, 2> i2;
    span_size = 1 << lvl;
    //std::cout<<"Spansize is :"<< span_size<<"\n";

    // left bottom corner box
    i1.setZero();
    i2[0] = bas->getBases()[0]->component(0).knots().unique().size()-1;
    i2[1] = bas->getBases()[0]->component(1).knots().unique().size()-1;

    i1[0] = i1[0] << (lvl); // get indices in lvl
    i1[1] = i1[1] << (lvl);
    i2[0] = i2[0] << (lvl);
    i2[1] = i2[1] << (lvl);

    boxes.push_back(i1);//insert the root of the tree into the boxes
    boxes.push_back(i2);
    srand((unsigned)time(NULL));//seed the random alg.
    for(int i = 0; i < lvl; i++){//insert lvl levels
        //std::cout<<"level: "<<i+1<< "\n";
        boxes_in_lvl = (rand()%max_nb);
        if(boxes_in_lvl == 0){
            boxes_in_lvl++;
        }
        for(int j = 0; j < boxes_in_lvl; j++){//insert number of boxes
            bi = rand() % (boxes.size()/2);//pick a box from prewious level
            if(aligned){
                //test if the box is a box not a point- insert only boxes
                i1[0] = ( ( (rand() % (boxes[2*bi+1][0]-boxes[2*bi][0]) ) + boxes[2*bi][0]) /span_size) << i;// left bottom corner box
                i1[1] = ( ( (rand() % (boxes[2*bi+1][1]-boxes[2*bi][1]) ) + boxes[2*bi][1]) /span_size) << i;
                i2[0] = ( ( (rand() % (boxes[2*bi+1][0]-i1[0]) +1) + i1[0]) /span_size) << i;// right top corner of the box
                i2[1] = ( ( (rand() % (boxes[2*bi+1][1]-i1[1]) +1) + i1[1]) /span_size) << i;

                i2[0] = math::max( i2[0], i1[0] + (bas->degree(1)+2)*(1 << lvl) );
                i2[1] = math::max( i2[1], i1[1] + (bas->degree(1)+2)*(1 << lvl) );
                
 
                //std::cout<<"\naligned box inserted ["<< i1[0]<<" , "<< i1[1]<<"] ["<< i2[0]<<" , "<< i2[1]<<"]"<< " to level "<<i+1 <<"\n";
                boxes_new.push_back(i1);
                boxes_new.push_back(i2);
                q.push_back(i+1);
                i1[0] = i1[0] >> (lvl - (i+1)); // get indices in level i+1
                i1[1] = i1[1] >> (lvl - (i+1));
                i2[0] = i2[0] >> (lvl - (i+1));
                i2[1] = i2[1] >> (lvl - (i+1));
                q.push_back(i1[0]);
                q.push_back(i1[1]);
                q.push_back(i2[0]);
                q.push_back(i2[1]);
                //bas->insert_box(i1,i2,i+1);
            }else{
                i1[0] = (rand() % (boxes[2*bi+1][0]-boxes[2*bi][0]) ) + boxes[2*bi][0];
                i1[1] = (rand() % (boxes[2*bi+1][1]-boxes[2*bi][1]) ) + boxes[2*bi][1];
                i2[0] = (rand() % (boxes[2*bi+1][0]-i1[0]) +1) + i1[0];
                i2[1] = (rand() % (boxes[2*bi+1][1]-i1[1]) +1) + i1[1];
                if( (i1[0]==i2[0]) || (i1[1]==i2[1]) ){
                    i2[0] += span_size;
                    i2[1] += span_size;
                }
                //std::cout<<"\nbox inserted ["<< i1[0]<<" , "<< i1[1]<<"] ["<< i2[0]<<" , "<< i2[1]<<"]"<< " to level "<<i+1 <<"\n";
                boxes_new.push_back(i1);
                boxes_new.push_back(i2);
                q.push_back(i+1);
                i1[0] = i1[0] >> (lvl - (i+1)); // get indices in level i+1
                i1[1] = i1[1] >> (lvl - (i+1));
                i2[0] = i2[0] >> (lvl - (i+1));
                i2[1] = i2[1] >> (lvl - (i+1));
                q.push_back(i1[0]);
                q.push_back(i1[1]);
                q.push_back(i2[0]);
                q.push_back(i2[1]);
                //bas->insert_box(i1,i2,i+1);
            }
        }
        span_size = span_size/2;
        //cout<<"Spansize is :"<< span_size<<endl;
        boxes.clear();
        for(unsigned int k = 0; k < boxes_new.size();k++){
            boxes.push_back(boxes_new[k]);
        }
        boxes_new.clear();
    }
    return q;
}
