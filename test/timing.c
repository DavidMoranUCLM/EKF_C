#include "EKF.h"
#include "esp_log.h"
#include "tictoc.h"
#include "esp_heap_trace.h"
#include "unity.h"

#define NUM_RECORDS 512
#define TIME_STEP_S 0.01f

static heap_trace_record_t trace_record[NUM_RECORDS];
const char *TAG = "EKF_C";
EKF_ctx_t ekf = {0};
measures_t meas = {0};

tictoc_t *initTic, *stepTic;



TEST_CASE("LEAKS", "[EKF_C]") {
  ESP_ERROR_CHECK(heap_trace_init_standalone(trace_record, NUM_RECORDS));
  ESP_ERROR_CHECK(heap_trace_start(HEAP_TRACE_LEAKS));

  ekfInit(&ekf, &meas);
  for(uint16_t i=0;i<100;i++){
    ekfStep(&ekf, &meas, i*TIME_STEP_S);
  }
  
  ekfDeinit(&ekf);

  heap_trace_stop();
  heap_trace_dump();

}



TEST_CASE("TIMING", "[EKF_C]") {
  ESP_LOGI(TAG,"Beginning test...");
  for (int i = 0; i < 100; i++) {
    //ESP_LOGI(TAG,"Iteration %d", i);
    tic(initTic);
    ekfInit(&ekf, &meas);
    toc(initTic);
    ekfDeinit(&ekf);
    ESP_LOGI(TAG, "Free heap size: %lu B", esp_get_free_heap_size());
  }
  tictoc_print(initTic);
}

void setUp(void) {
  meas.acc[0] = 0;
  meas.acc[1] = 2;
  meas.acc[2] = 9.5;
  meas.mag[0] = 0;
  meas.mag[1] = 1;
  meas.mag[2] = 0;
  initTic = tictoc_new("Init time");
  stepTic = tictoc_new("Step time");
  ESP_LOGI(TAG,"Timers Set");
}