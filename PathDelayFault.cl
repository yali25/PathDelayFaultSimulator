__kernel void
logicsimulationv1(ulong length, __global uint *top_order,
                  __global uint *reverse_top_order, __global uint *gate_type,
                  __global uint2 *edges, __global uint2 *out_edges,
                  __global uint *p1, __global uint *p2,
                  __global uint *robust_result1, __global uint *robust_result2,
                  __global uint *is_robust_or_non_robust_result1,
                  __global uint *is_robust_or_non_robust_result2,
                  __global uint *is_robust_path_to_po,
                  __global uint *belongs_to_robust_path_final_result1,
                  __global uint *belongs_to_robust_path_final_result2) {
  const unsigned int output_offset = length * get_global_id(0);
  // const ulong size = length;
  unsigned int result1 = 0;
  unsigned int result2 = 0;
  unsigned int pattern1_1;
  unsigned int pattern1_2;
  unsigned int pattern2_1;
  unsigned int pattern2_2;
  unsigned int isRobust_1;
  unsigned int isRobust_2;
  unsigned int isRobustOrNonRobust_l;
  unsigned int isRobustOrNonRobust_2;
  for (int i = 0; i < SIZE; i++) {
    const unsigned int c_node = i;
    result1 = 0;
    result2 = 0;
    const unsigned int g = gate_type[c_node];
    if (g) // All other gates, if node is an input the else is executed, the
           // gate type can also be used to determine the in_degree
    {

      const uint2 e = edges[c_node];
      pattern1_1 = p1[output_offset + e.x];
      pattern1_2 = p1[output_offset + e.y];
      pattern2_1 = p2[output_offset + e.x];
      pattern2_2 = p2[output_offset + e.y];

      if (g == 1) {

        result1 = ~(pattern1_1 & pattern1_2);
        result2 = ~(pattern2_1 & pattern2_2);
      } else if (g == 2) {

        result1 = ~(pattern1_1 | pattern1_2);
        result2 = ~(pattern2_1 | pattern2_2);
      } else if (g == 4) {

        result1 = pattern1_1 & pattern1_2;
        result2 = pattern2_1 & pattern2_2;
      } else if (g == 5) {

        result1 = pattern1_1 | pattern1_2;
        result2 = pattern2_1 | pattern2_2;
      } else if (g == 6) {
        result1 = ~pattern1_1;
        result2 = ~pattern2_1;
      }

      if (g == 1 || g == 4) // AND NAND case
      {
        // isRobustAND NAND
        isRobust_1 = ((~pattern1_1) & pattern1_2 & pattern2_2) |
                     (pattern1_1 & (~pattern1_2) & pattern2_1 & pattern2_2);
        isRobust_2 = ((~pattern2_1) & pattern2_2 & pattern1_2) |
                     (pattern2_1 & (~pattern2_2) & pattern1_1 & pattern1_2);

        isRobustOrNonRobust_l = (pattern1_1 ^ pattern1_2) & pattern2_2;
        isRobustOrNonRobust_2 = (pattern2_1 ^ pattern2_2) & pattern1_2;
      } else if (g == 5 || g == 2) // Or Nor
      {
        isRobust_1 =
            (pattern1_1 & (~pattern1_2) & (~pattern2_2)) |
            ((~pattern1_1) & pattern1_2 & (~pattern2_1) & (~pattern2_2));
        isRobust_2 =
            (pattern2_1 & (~pattern2_2) & (~pattern1_2)) |
            ((~pattern2_1) & pattern2_2 & (~pattern1_1) & (~pattern1_2));

        isRobustOrNonRobust_l = (pattern1_1 ^ pattern1_2) & (~pattern2_2);
        isRobustOrNonRobust_2 = (pattern2_1 ^ pattern2_2) & (~pattern1_2);
      } else if (g == 6) {
        isRobust_1 = pattern1_1 ^ pattern2_1;
        isRobust_2 = pattern1_1 ^ pattern2_1;

        isRobustOrNonRobust_l = pattern1_1 ^ pattern2_1;
        isRobustOrNonRobust_2 = pattern1_1 ^ pattern2_1;
      }
      p1[output_offset + c_node] = result1;
      p2[output_offset + c_node] = result2;

      robust_result1[output_offset + c_node] = isRobust_1;
      robust_result2[output_offset + c_node] = isRobust_2;
      is_robust_or_non_robust_result1[output_offset + c_node] =
          isRobustOrNonRobust_l;
      is_robust_or_non_robust_result2[output_offset + c_node] =
          isRobustOrNonRobust_2;
    }
  }

  /* traverse over the graph in reverse topological order */
  //#pragma unroll 32
  for (int i = SIZE - 1; i >= 0; i--) {
    const unsigned int c_node = i;
    const uint2 e = out_edges[c_node];
    if ((e.x != 0) && (e.y != 0)) // Node is not an output
    // if (out_edges[c_node])
    {
      unsigned int result = 0;

      unsigned int r1 = 0;
      unsigned int path_po = 0;
      r1 = (edges[out_edges[c_node].x].x == c_node)
               ? robust_result1[output_offset + out_edges[c_node].x]
               : robust_result2[output_offset + out_edges[c_node].x];
      path_po = is_robust_path_to_po[output_offset + out_edges[c_node].x];
      result |= (r1 & path_po);

      r1 = (edges[out_edges[c_node].y].x == c_node)
               ? robust_result1[output_offset + out_edges[c_node].y]
               : robust_result2[output_offset + out_edges[c_node].y];
      path_po = is_robust_path_to_po[output_offset + out_edges[c_node].y];
      result |= (r1 & path_po);
      is_robust_path_to_po[output_offset + c_node] = result;
    } else {
      is_robust_path_to_po[output_offset + c_node] =
          ~0; // set the bits of all outputs to 1
    }         // Node is an output because no out edges
              // case of output set all bits to "1"
  }

  /* Last Interation over the circuit in topological order*/
  //#pragma unroll 32
  for (int i = 0; i < SIZE; i++) {
    const unsigned int c_node = i;
    const unsigned int g = gate_type[c_node];
    if (g) {
      /* Check all input edges */
      unsigned int temp_final_result = 0;
      unsigned int btrofr = 0;

      btrofr = (out_edges[edges[c_node].x].x == c_node)
                   ? belongs_to_robust_path_final_result1[output_offset +
                                                          edges[c_node].x]
                   : belongs_to_robust_path_final_result2[output_offset +
                                                          edges[c_node].x];
      temp_final_result |=
          is_robust_path_to_po[output_offset + edges[c_node].x] &
          robust_result1[output_offset + c_node] & btrofr;

      btrofr = (out_edges[edges[c_node].y].x == c_node)
                   ? belongs_to_robust_path_final_result1[output_offset +
                                                          edges[c_node].y]
                   : belongs_to_robust_path_final_result2[output_offset +
                                                          edges[c_node].y];
      temp_final_result |=
          is_robust_path_to_po[output_offset + edges[c_node].y] &
          robust_result2[output_offset + c_node] & btrofr;

      /* Check all output edges */
      if ((out_edges[c_node].x != 0) && (out_edges[c_node].y != 0)) {
        result1 = 0;
        result2 = 0;

        unsigned int t = 0;

        t = (edges[out_edges[c_node].x].x == c_node)
                ? robust_result1[output_offset + out_edges[c_node].x]
                : robust_result2[output_offset + out_edges[c_node].x];
        result1 = t & temp_final_result &
                  is_robust_path_to_po[output_offset + out_edges[c_node].x];

        t = (edges[out_edges[c_node].y].x == c_node)
                ? robust_result1[output_offset + out_edges[c_node].y]
                : robust_result2[output_offset + out_edges[c_node].y];
        result2 = t & temp_final_result &
                  is_robust_path_to_po[output_offset + out_edges[c_node].y];

        belongs_to_robust_path_final_result1[output_offset + c_node] = result1;
        belongs_to_robust_path_final_result2[output_offset + c_node] = result2;

      } else {
        belongs_to_robust_path_final_result1[output_offset + c_node] =
            temp_final_result;
        belongs_to_robust_path_final_result2[output_offset + c_node] =
            temp_final_result;
      }
    } else {
      belongs_to_robust_path_final_result1[output_offset + c_node] =
          is_robust_path_to_po[output_offset + c_node];
      belongs_to_robust_path_final_result2[output_offset + c_node] =
          is_robust_path_to_po[output_offset + c_node];
    }
  }
}

__kernel void
logicsimulationv2(__global uint *gate_type, __global uint2 *edges,
                  __global uint2 *out_edges, __global uint *p1,
                  __global uint *p2, __global uint *robust_result1,
                  __global uint *robust_result2,
                  __global uint *is_robust_or_non_robust_result1,
                  __global uint *is_robust_or_non_robust_result2,
                  __global uint *is_robust_path_to_po,
                  __global uint *belongs_to_robust_path_final_result1,
                  __global uint *belongs_to_robust_path_final_result2) {
  const unsigned int output_offset = SIZE * get_global_id(0);
  unsigned int result1 = 0;
  unsigned int result2 = 0;
  unsigned int pattern1_1;
  unsigned int pattern1_2;
  unsigned int pattern2_1;
  unsigned int pattern2_2;
  unsigned int isRobust_1;
  unsigned int isRobust_2;
  unsigned int isRobustOrNonRobust_l;
  unsigned int isRobustOrNonRobust_2;
  for (int i = 0; i < SIZE; i++) {
    const unsigned int c_node = i;
    result1 = 0;
    result2 = 0;
    const unsigned int g = gate_type[c_node];
    if (g) // All other gates, if node is an input the else is executed, the
           // gate type can also be used to determine the in_degree
    {

      const uint2 e = edges[c_node];
      pattern1_1 = p1[output_offset + e.x];
      pattern1_2 = p1[output_offset + e.y];
      pattern2_1 = p2[output_offset + e.x];
      pattern2_2 = p2[output_offset + e.y];

      if (g == 1) {

        result1 = ~(pattern1_1 & pattern1_2);
        result2 = ~(pattern2_1 & pattern2_2);
      } else if (g == 2) {

        result1 = ~(pattern1_1 | pattern1_2);
        result2 = ~(pattern2_1 | pattern2_2);
      } else if (g == 4) {

        result1 = pattern1_1 & pattern1_2;
        result2 = pattern2_1 & pattern2_2;
      } else if (g == 5) {

        result1 = pattern1_1 | pattern1_2;
        result2 = pattern2_1 | pattern2_2;
      } else if (g == 6) {
        result1 = ~pattern1_1;
        result2 = ~pattern2_1;
      }

      if (g == 1 || g == 4) // AND NAND case
      {
        // isRobustAND NAND
        isRobust_1 = ((~pattern1_1) & pattern1_2 & pattern2_2) |
                     (pattern1_1 & (~pattern1_2) & pattern2_1 & pattern2_2);
        isRobust_2 = ((~pattern2_1) & pattern2_2 & pattern1_2) |
                     (pattern2_1 & (~pattern2_2) & pattern1_1 & pattern1_2);

        isRobustOrNonRobust_l = (pattern1_1 ^ pattern1_2) & pattern2_2;
        isRobustOrNonRobust_2 = (pattern2_1 ^ pattern2_2) & pattern1_2;
      } else if (g == 5 || g == 2) // Or Nor
      {
        isRobust_1 =
            (pattern1_1 & (~pattern1_2) & (~pattern2_2)) |
            ((~pattern1_1) & pattern1_2 & (~pattern2_1) & (~pattern2_2));
        isRobust_2 =
            (pattern2_1 & (~pattern2_2) & (~pattern1_2)) |
            ((~pattern2_1) & pattern2_2 & (~pattern1_1) & (~pattern1_2));

        isRobustOrNonRobust_l = (pattern1_1 ^ pattern1_2) & (~pattern2_2);
        isRobustOrNonRobust_2 = (pattern2_1 ^ pattern2_2) & (~pattern1_2);
      } else if (g == 6) {
        isRobust_1 = pattern1_1 ^ pattern2_1;
        isRobust_2 = pattern1_1 ^ pattern2_1;

        isRobustOrNonRobust_l = pattern1_1 ^ pattern2_1;
        isRobustOrNonRobust_2 = pattern1_1 ^ pattern2_1;
      }
      const uint o = output_offset + c_node;
      p1[o] = result1;
      p2[o] = result2;

      robust_result1[o] = isRobust_1;
      robust_result2[o] = isRobust_2;
      is_robust_or_non_robust_result1[o] = isRobustOrNonRobust_l;
      is_robust_or_non_robust_result2[o] = isRobustOrNonRobust_2;
    }
  }

  /* traverse over the graph in reverse topological order */
  //#pragma unroll 32
  for (int i = SIZE - 1; i >= 0; i--) {
    const unsigned int c_node = i;
    const uint2 e = out_edges[c_node];
    if ((e.x != 0) && (e.y != 0)) // Node is not an output
    // if (out_edges[c_node])
    {
      unsigned int result = 0;

      unsigned int r1 = 0;
      unsigned int path_po = 0;
      r1 = (edges[out_edges[c_node].x].x == c_node)
               ? robust_result1[output_offset + out_edges[c_node].x]
               : robust_result2[output_offset + out_edges[c_node].x];
      path_po = is_robust_path_to_po[output_offset + out_edges[c_node].x];
      result |= (r1 & path_po);

      r1 = (edges[out_edges[c_node].y].x == c_node)
               ? robust_result1[output_offset + out_edges[c_node].y]
               : robust_result2[output_offset + out_edges[c_node].y];
      path_po = is_robust_path_to_po[output_offset + out_edges[c_node].y];
      result |= (r1 & path_po);
      is_robust_path_to_po[output_offset + c_node] = result;
    } else {
      is_robust_path_to_po[output_offset + c_node] =
          ~0; // set the bits of all outputs to 1
    }         // Node is an output because no out edges
              // case of output set all bits to "1"
  }

  /* Last Interation over the circuit in topological order*/
  //#pragma unroll 32
  for (int i = 0; i < SIZE; i++) {
    const unsigned int c_node = i;
    const unsigned int g = gate_type[c_node];
    if (g) {
      /* Check all input edges */
      unsigned int temp_final_result = 0;
      unsigned int btrofr = 0;

      btrofr = (out_edges[edges[c_node].x].x == c_node)
                   ? belongs_to_robust_path_final_result1[output_offset +
                                                          edges[c_node].x]
                   : belongs_to_robust_path_final_result2[output_offset +
                                                          edges[c_node].x];
      temp_final_result |=
          is_robust_path_to_po[output_offset + edges[c_node].x] &
          robust_result1[output_offset + c_node] & btrofr;

      btrofr = (out_edges[edges[c_node].y].x == c_node)
                   ? belongs_to_robust_path_final_result1[output_offset +
                                                          edges[c_node].y]
                   : belongs_to_robust_path_final_result2[output_offset +
                                                          edges[c_node].y];
      temp_final_result |=
          is_robust_path_to_po[output_offset + edges[c_node].y] &
          robust_result2[output_offset + c_node] & btrofr;

      /* Check all output edges */
      if ((out_edges[c_node].x != 0) && (out_edges[c_node].y != 0)) {
        result1 = 0;
        result2 = 0;

        unsigned int t = 0;

        t = (edges[out_edges[c_node].x].x == c_node)
                ? robust_result1[output_offset + out_edges[c_node].x]
                : robust_result2[output_offset + out_edges[c_node].x];
        result1 = t & temp_final_result &
                  is_robust_path_to_po[output_offset + out_edges[c_node].x];

        t = (edges[out_edges[c_node].y].x == c_node)
                ? robust_result1[output_offset + out_edges[c_node].y]
                : robust_result2[output_offset + out_edges[c_node].y];
        result2 = t & temp_final_result &
                  is_robust_path_to_po[output_offset + out_edges[c_node].y];

        belongs_to_robust_path_final_result1[output_offset + c_node] = result1;
        belongs_to_robust_path_final_result2[output_offset + c_node] = result2;

      } else {
        belongs_to_robust_path_final_result1[output_offset + c_node] =
            temp_final_result;
        belongs_to_robust_path_final_result2[output_offset + c_node] =
            temp_final_result;
      }
    } else {
      belongs_to_robust_path_final_result1[output_offset + c_node] =
          is_robust_path_to_po[output_offset + c_node];
      belongs_to_robust_path_final_result2[output_offset + c_node] =
          is_robust_path_to_po[output_offset + c_node];
    }
  }
}
